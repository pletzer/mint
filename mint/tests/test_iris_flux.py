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
DEG2RAD = np.pi / 180.

def _set_extensive_field_from_streamfct(cube):
    """
    Set extensive field values from streamfunction cos(y)*cos(x)
    :param cube: the cube for which we will fill in the values
    """
    x = cube.mesh.node_coords.node_x.points
    y = cube.mesh.node_coords.node_y.points

    # make indexing is zero-based
    e2n = cube.mesh.edge_node_connectivity.indices_by_location() - \
          cube.mesh.edge_node_connectivity.start_index

    num_edges = e2n.shape[0]
    for edge in range(num_edges):
        n0, n1 = e2n[edge, :]
        # the end points of the edge
        x0, x1 = x[n0], x[n1]
        y0, y1 = y[n0], y[n1]
        s0 = np.cos(y0*DEG2RAD)*np.cos(x0*DEG2RAD)
        s1 = np.cos(y1*DEG2RAD)*np.cos(x1*DEG2RAD)
        cube.data[edge] = s1 - s0

def _set_vector_field_from_potentialfct(u_cube, v_cube):
    """
    Set vector field values from potential cos(y)*cos(x)
    :param u_cube: the x-component cube for which we will fill in the values
    :param v_cube: the y-component cube for which we will fill in the values
    """
    xe, ye, ze = mint.computeEdgeXYZ(u_cube.mesh, radius=1.0)
    lone, late = mint.computeLonLatFromXYZ(xe, ye, ze)
    num_edges = len(lone)
    other_dims = u_cube.shape[:-1]
    mai = mint.MultiArrayIter(other_dims)
    mai.begin()
    for _ in range(mai.getNumIters()):
        inds = tuple(mai.getIndices())
        slab = inds + (slice(0, num_edges),)
        u_cube.data[:] = - np.sin(lone*DEG2RAD)
        v_cube.data[:] = - np.sin(late*DEG2RAD) * np.sin(lone*DEG2RAD)


def _set_vector_field_from_streamfct(u_cube, v_cube):
    """
    Set vector field values from streamfunction cos(y)*cos(x)
    :param u_cube: the x-component cube for which we will fill in the values
    :param v_cube: the y-component cube for which we will fill in the values
    """
    xe, ye, ze = mint.computeEdgeXYZ(u_cube.mesh, radius=1.0)
    lone, late = mint.computeLonLatFromXYZ(xe, ye, ze)
    num_edges = len(lone)
    other_dims = u_cube.shape[:-1]
    mai = mint.MultiArrayIter(other_dims)
    mai.begin()
    for _ in range(mai.getNumIters()):
        inds = tuple(mai.getIndices())
        slab = inds + (slice(0, num_edges),)
        u_cube.data[slab] = - np.sin(late*DEG2RAD) * np.cos(lone*DEG2RAD)
        v_cube.data[slab] = np.sin(lone*DEG2RAD)
        mai.next()



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


    shape = tuple([])
    extra_coords = []
    if time is not None:
        shape += (time,)
        time_coord = DimCoord(np.arange(time), standard_name="time", units="days since 1970-01-01")
        time_coord.guess_bounds()
        extra_coords.append((time_coord, len(extra_coords)))
    if height is not None:
        shape += (height,)
        height_coord = DimCoord(np.arange(height), standard_name="height", units="meters")
        height_coord.guess_bounds()
        extra_coords.append((height_coord, len(extra_coords)))

    shape += mesh_coord_x.points.shape

    data = np.zeros(shape)
    cube = Cube(data)
    for coord, dim in extra_coords:
        cube.add_dim_coord(coord, dim)
    cube.add_aux_coord(mesh_coord_x, len(extra_coords))
    cube.add_aux_coord(mesh_coord_y, len(extra_coords))
    # cube has a mesh (cube.mesh)
    return cube


def test_cs_zt():    

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind_zt.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)

    xy = [(0., 0.), (360., 0.), (160., 70.), (0., 0.)]
    flx_calc = mint.IrisMintFlux(src_u.mesh, src_flags=src_flags, tgt_line=xy)
    flux = flx_calc.evaluate_from_vector(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    # check
    assert np.all(np.fabs(flux - 0.0) < 1.e-10)




