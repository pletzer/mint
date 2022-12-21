import iris
from iris.coords import AuxCoord, DimCoord
from iris.cube import Cube
from iris.experimental.ugrid import Connectivity, Mesh, PARSE_UGRID_ON_LOAD
import numpy as np
from numpy import ma
from pathlib import Path

from mint.iris_regrid import MINTScheme
import mint



DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

def _set_extensive_field_from_streamfct(cube):
    """
    Set extensive field values from streamfunction cos(y)*cos(x)
    :param cube: the cube for which we will fill in the values
    """
    x = cube.mesh.node_coords.node_x.points
    y = cube.mesh.node_coords.node_y.points

    e2n = cube.mesh.edge_node_connectivity.indices_by_location()
    # make sure the edge to node connectivity is zero based
    e2n -= cube.mesh.edge_node_connectivity.start_index

    num_edges = e2n.shape[0]
    deg2rad = np.pi/180.
    for edge in range(num_edges):
        n0, n1 = e2n[edge, :]
        # the end points of the edge
        x0, x1 = x[n0], x[n1]
        y0, y1 = y[n0], y[n1]
        s0 = np.cos(y0*deg2rad)*np.cos(x0*deg2rad)
        s1 = np.cos(y1*deg2rad)*np.cos(x1*deg2rad)
        cube.data[edge] = s1 - s0

def _set_vector_field_from_potentialfct(u_cube, v_cube):
    """
    Set vector field values from potential cos(y)*cos(x)
    :param u_cube: the x-component cube for which we will fill in the values
    :param v_cube: the y-component cube for which we will fill in the values
    """
    x = u_cube.mesh.node_coords.node_x.points
    y = u_cube.mesh.node_coords.node_y.points

    e2n = u_cube.mesh.edge_node_connectivity.indices_by_location()
    # make sure the edge to node connectivity is zero based
    e2n -= u_cube.mesh.edge_node_connectivity.start_index

    num_edges = e2n.shape[0]
    deg2rad = np.pi/180.
    for edge in range(num_edges):
        n0, n1 = e2n[edge, :]
        # the end points of the edge
        x0, x1 = x[n0], x[n1]
        y0, y1 = y[n0], y[n1]
        # mid point on the edge
        xm = 0.5*(x0 + x1)
        ym = 0.5*(y0 + y1)
        u_cube.data[edge] = - np.sin(xm*deg2rad)
        v_cube.data[edge] = - np.sin(ym*deg2rad) * np.sin(xm*deg2rad)

def _set_vector_field_from_streamfct(u_cube, v_cube):
    """
    Set vector field values from streamfunction cos(y)*cos(x)
    :param u_cube: the x-component cube for which we will fill in the values
    :param v_cube: the y-component cube for which we will fill in the values
    """
    x = u_cube.mesh.node_coords.node_x.points
    y = u_cube.mesh.node_coords.node_y.points

    e2n = u_cube.mesh.edge_node_connectivity.indices_by_location()
    # make sure the edge to node connectivity is zero based
    e2n -= u_cube.mesh.edge_node_connectivity.start_index

    num_edges = e2n.shape[0]
    deg2rad = np.pi/180.
    for edge in range(num_edges):
        n0, n1 = e2n[edge, :]
        # the end points of the edge
        x0, x1 = x[n0], x[n1]
        y0, y1 = y[n0], y[n1]
        # mid point on the edge
        xm = 0.5*(x0 + x1)
        ym = 0.5*(y0 + y1)
        u_cube.data[edge] = - np.sin(ym*deg2rad) * np.cos(xm*deg2rad)
        v_cube.data[edge] = np.sin(xm*deg2rad)


def _u_v_cubes_from_ugrid_file(filename, 
                               u_std_name: str="eastward_wind_at_cell_faces",
                               v_std_name: str="northward_wind_at_cell_faces"):
    """
    Get u, v components cubes from a Ugrid file
    :param filename: netCDF file name
    :param u_std_name: standard name for the zonal component of the vector field
    :param v_std_name: standard name for the meridional component of the vector field
    :returns (u_cube, v_cube)
    """
    u_std_name = "eastward_wind_at_cell_faces"
    v_std_name = "eastward_wind_at_cell_faces"
    with PARSE_UGRID_ON_LOAD.context():
        u_cube = iris.load_cube(filename, u_std_name)
        v_cube = iris.load_cube(filename, v_std_name)
    return (u_cube, v_cube)


def _gridlike_mesh(n_lons, n_lats):
    """
    Generate a global mesh with geometry similar to a rectilinear grid.
    The resulting mesh will have n_lons cells spanning its longitudes and
    n_lats cells spanning its latitudes for a total of (n_lons * n_lats) cells.
    """
    n_lons1 = n_lons + 1
    n_lats1 = n_lats + 1

    # Latitude and longitude values are set.
    lat_values = np.linspace(-90, 90, n_lats1)
    lon_values = np.linspace(-180, 180, n_lons1)
    lon_array, lat_array = np.meshgrid(lon_values, lat_values)

    lon_array = lon_array.flatten()
    lat_array = lat_array.flatten()
    n_points = lon_array.shape[0]

    # Build the face to node connectivity.
    n_cells = n_lons * n_lats
    face_lons = np.empty((n_cells,), np.float64)
    face_lats = np.empty((n_cells,), np.float64)

    fnc_array = np.empty((n_cells, 4), np.int32)
    cell_index = 0
    for j0 in range(n_lats):
        j1 = j0 + 1
        for i0 in range(n_lons):
            i1 = i0 + 1
            node_inds = (i0 + n_lons1*j0, i1 + n_lons1*j0, i1 + n_lons1*j1, i0 + n_lons1*j1)
            fnc_array[cell_index, :] = node_inds
            face_lons[cell_index] = lon_array.take(node_inds).mean()
            face_lats[cell_index] = lat_array.take(node_inds).mean()
            cell_index += 1

    # Build the edge to node connectivity.
    n_edges = n_lons*(n_lats + 1) + (n_lons + 1)*n_lats
    edge_lons = np.empty((n_edges,), np.float64)
    edge_lats = np.empty((n_edges,), np.float64)
    enc_array = np.empty((n_edges, 2), np.int32)
    # horizontal edges
    edge_index = 0
    for j in range(n_lats1):
        for i0 in range(n_lons):
            i1 = i0 + 1
            node_inds = (i0 + n_lons1*j, i1 + n_lons1*j)
            enc_array[edge_index] = node_inds
            edge_lons[edge_index] = lon_array.take(node_inds).mean()
            edge_lats[edge_index] = lat_array.take(node_inds).mean()
            edge_index += 1
    for j0 in range(n_lats):
        j1 = j0 + 1
        for i in range(n_lons1):
            node_inds = (i + n_lons1*j0, i + n_lons1*j1)
            enc_array[edge_index] = node_inds
            edge_lons[edge_index] = lon_array.take(node_inds).mean()
            edge_lats[edge_index] = lat_array.take(node_inds).mean()
            edge_index += 1
    
    # Translate the mesh information into iris objects.
    fnc = Connectivity(
        fnc_array,
        cf_role="face_node_connectivity",
        start_index=0,
    )
    enc = Connectivity(
        enc_array,
        cf_role="edge_node_connectivity",
        start_index=0,
    )
    lons = AuxCoord(lon_array, standard_name="longitude")
    lats = AuxCoord(lat_array, standard_name="latitude")

    # In order to add a mesh to a cube, face locations must be added.
    face_lon_coord = AuxCoord(face_lons, standard_name="longitude")
    face_lat_coord = AuxCoord(face_lats, standard_name="latitude")

    # Add edge coords.
    edge_lon_coord = AuxCoord(edge_lons, standard_name="longitude")
    edge_lat_coord = AuxCoord(edge_lats, standard_name="latitude")

    mesh = Mesh(2, ((lons, "x"), (lats, "y")), [fnc, enc])

    mesh.add_coords(
        face_x=face_lon_coord,
        face_y=face_lat_coord,
        edge_x=edge_lon_coord,
        edge_y=edge_lat_coord,
    )
    mesh.long_name = "example mesh"

    mint.saveMeshVTK(mesh, 'example_mesh.vtk')

    return mesh


def _gridlike_mesh_cube(n_lons, n_lats, location="edge", time=None, height=None):
    """
    Create a grid-like mesh cube
    :param n_lons: number of lon cells
    :param n_lats: number of lat cells
    :param location: edge mesh?
    :param time: number of time intervals (int)
    :param height: number of elevation intervals (int)
    """
    mesh = _gridlike_mesh(n_lons, n_lats)
    mesh_coord_x, mesh_coord_y = mesh.to_MeshCoords(location)
    shape = mesh_coord_x.points.shape
    extra_coords = []
    if time is not None:
        shape += (time,)
        time_coord = DimCoord(np.arange(time), standard_name="time", units="days since 1970-01-01")
        time_coord.guess_bounds()
        extra_coords.append((time_coord, len(extra_coords) + 1))
    if height is not None:
        shape += (height,)
        height_coord = DimCoord(np.arange(height), standard_name="height", units="meters")
        height_coord.guess_bounds()
        extra_coords.append((height_coord, len(extra_coords) + 1))
    data = np.zeros(shape)
    cube = Cube(data)
    cube.add_aux_coord(mesh_coord_x, 0)
    cube.add_aux_coord(mesh_coord_y, 0)
    for coord, dim in extra_coords:
        cube.add_dim_coord(coord, dim)
    # cube has a mesh (cube.mesh)
    return cube


def test_read_ugrid_file():
    u_cube, v_cube = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs8_wind.nc')


def test_cubedsphere8_to_cubedsphere2():

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))
    tgt_u, tgt_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs2_wind.nc'))

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)
    # grid options for a cubed-sphere grid
    tgt_flags = (1, 1, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags)

    _set_vector_field_from_streamfct(src_u, src_v)
    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    # Check.
    error = 0.5*(np.fabs(result_u.data - tgt_u.data).mean() + \
                 np.fabs(result_v.data - tgt_v.data).mean())
    print(f'test_cubedsphere8_to_cubedsphere2 = {error}')
    assert error < 0.29


# turned test off as we have not yet implemented W1
def test_cubedsphere8_to_cubedsphere8_w1():

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))
    tgt_u, tgt_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)
    # grid options for a cubed-sphere grid
    tgt_flags = (1, 1, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags)

    _set_vector_field_from_potentialfct(src_u, src_v)
    _set_vector_field_from_potentialfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W1)

    # Check.
    error = 0.5*(np.fabs(result_u.data - tgt_u.data).mean() + \
                 np.fabs(result_v.data - tgt_v.data).mean())
    print(f'test_cubedsphere8_to_cubedsphere8_w1 = {error}')
    assert error < 0.06

def test_cubedsphere8_to_cubedsphere8_w2():

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))
    tgt_u, tgt_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)
    # grid options for a cubed-sphere grid
    tgt_flags = (1, 1, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags)

    _set_vector_field_from_streamfct(src_u, src_v)
    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    mint.saveVectorFieldVTK(result_u, result_v, 'result_vectors.vtk')

    # Check.
    error = 0.5*(np.fabs(result_u.data - tgt_u.data).mean() + \
                 np.fabs(result_v.data - tgt_v.data).mean())
    print(f'test_cubedsphere8_to_cubedsphere8_w2 = {error}')
    assert error < 0.06


def test_lonlat_to_cubedsphere():

    src_u = _gridlike_mesh_cube(9, 5)
    src_v = _gridlike_mesh_cube(9, 5)

    tgt_u, tgt_v = _u_v_cubes_from_ugrid_file(DATA_DIR / Path('cs8_wind.nc'))

    # grid options for a regular lat-lon grid
    src_flags = (0, 0, 1)
    # grid options for a cubed-sphere grid
    tgt_flags = (1, 1, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags)

    _set_vector_field_from_streamfct(src_u, src_v)
    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    # Check.
    error = 0.5*(np.fabs(result_u.data - tgt_u.data).mean() + \
                 np.fabs(result_v.data - tgt_v.data).mean())
    print(f'test_lonlat_to_cubedsphere = {error}')

    assert error < 0.10


def test_cube_mesh():

    cube = _gridlike_mesh_cube(4, 5)
    assert hasattr(cube, 'shape')
    assert hasattr(cube, 'data')
    assert hasattr(cube, 'mesh')


def test_mesh_to_mesh_basic():

    src = _gridlike_mesh_cube(4, 5)
    tgt = _gridlike_mesh_cube(6, 3)
    src_mesh = src.mesh
    tgt_mesh = tgt.mesh

    # compute the regridding weights
    rg = mint.IrisMintRegridder(src_mesh, tgt_mesh,
                                src_flags=(0, 0, 1), tgt_flags=(0, 0, 1))

    # extensive field regridding
    out_cube = rg.regrid_extensive_cube(src)


def test_streamfunction_extensive_field():

    # Create source cubes on unstructured meshes.
    src_nx, src_ny = 20, 10
    src_nx1, src_ny1 = src_nx + 1, src_ny + 1
    src_num_edges = src_nx*src_ny1 + src_nx1*src_ny
    src = _gridlike_mesh_cube(src_nx, src_ny)

    src = _gridlike_mesh_cube(20, 10)
    tgt = _gridlike_mesh_cube(30, 20)
    
    rg = mint.IrisMintRegridder(src.mesh, tgt.mesh, \
                                src_flags=(0, 0, 1), tgt_flags=(0, 0, 1))

    # Regrid the extensive field from stream function.
    _set_extensive_field_from_streamfct(src)
    _set_extensive_field_from_streamfct(tgt)
    result = rg.regrid_extensive_cube(src)

    # Check the result.
    error = np.mean(np.fabs(result.data - tgt.data))
    print(f'extensive field regridding error = {error}')
    assert error < 0.007


def test_streamfunction_vector_field():

    # Regrid the vector field from stream function.
    src_nx, src_ny = 18, 14
    src_nx1, src_ny1 = src_nx + 1, src_ny + 1
    src_num_edges = src_nx*src_ny1 + src_nx1*src_ny
    src_u = _gridlike_mesh_cube(src_nx, src_ny)
    src_v = _gridlike_mesh_cube(src_nx, src_ny)
    _set_vector_field_from_streamfct(src_u, src_v)
    assert src_u.shape[-1] == src_v.shape[-1] == src_num_edges

    tgt_nx, tgt_ny = 14, 12
    tgt_nx1, tgt_ny1 = tgt_nx + 1, tgt_ny + 1
    tgt_num_edges = tgt_nx*tgt_ny1 + tgt_nx1*tgt_ny
    tgt_u = _gridlike_mesh_cube(tgt_nx, tgt_ny)
    tgt_v = _gridlike_mesh_cube(tgt_nx, tgt_ny)
    _set_vector_field_from_streamfct(tgt_u, tgt_v)
    assert tgt_u.shape[-1] == tgt_v.shape[-1] == tgt_num_edges

    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=(0, 0, 1), tgt_flags=(0, 0, 1))
    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)
    assert result_u.shape[-1] == result_v.shape[-1] == tgt_num_edges

    print(f'exact u = {tgt_u.data}')
    print(f'regridded u = {result_u.data}')

    # Check the result.
    tgt_x = result_u.mesh.node_coords.node_x.points
    tgt_y = result_u.mesh.node_coords.node_y.points
    tgt_num_cells = rg.tgt.get_grid().getNumberOfCells()
    error = 0.0
    eps = 1.e-8
    for icell in range(tgt_num_cells):
        for ie in range(mint.NUM_EDGES_PER_QUAD):
            # get the edge points
            point_id0, point_id1 = rg.tgt.get_grid().getNodeIds(icell, ie)
            x0, y0 = tgt_x[point_id0], tgt_y[point_id0]
            x1, y1 = tgt_x[point_id1], tgt_y[point_id1]
            edge_id, edge_sign = rg.tgt.get_grid().getEdgeId(icell, ie)
            if abs(abs(y0) - 90) < eps and abs(abs(y1) - 90) < eps:
                # the 2 nodes are at the pole, u, v are ill-defined thre so skip
                continue
            uval, vval = result_u.data[edge_id], result_v.data[edge_id]
            uxct, vxct = tgt_u.data[edge_id], tgt_v.data[edge_id]
            error += abs(uval - uxct)
            error += abs(vval - vxct)
            print(f'cell {icell} edge {ie} points {x0, y0} -> {x1, y1} u,v = {uval, vval} (exact {uxct, vxct})')
    error /= (2 * tgt_num_edges)
    print(f'error = {error}')
    assert error < 0.04

