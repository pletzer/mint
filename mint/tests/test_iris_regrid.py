from iris.coords import AuxCoord
from iris.cube import Cube
from iris.experimental.ugrid import Connectivity, Mesh
import numpy as np
from numpy import ma

from mint.iris_regrid import MINTScheme
import mint


def _set_extensive_field_data_from_streamfct(cube):
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
        s0 = np.sin(y0*deg2rad) + np.cos(y0*deg2rad)*np.cos(x0*deg2rad)
        s1 = np.sin(y1*deg2rad) + np.cos(y1*deg2rad)*np.cos(x1*deg2rad)
        cube.data[edge] = s1 - s0

def _gridlike_mesh2(n_lons, n_lats):
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

    fnc_array = np.empty((n_cells, 4), np.int)
    cell_index = 0
    for j0 in range(n_lats):
        j1 = j0 + 1
        for i0 in range(n_lons):
            i1 = i0 + 1
            node_inds = (i0 + n_lons*j0, i1 + n_lons*j0, i1 + n_lons*j1, i0 + n_lons*j1)
            fnc_array[cell_index, :] = node_inds
            face_lons[cell_index] = lon_array.take(node_inds).mean()
            face_lats[cell_index] = lat_array.take(node_inds).mean()
            cell_index += 1

    # Build the edge to node connectivity.
    n_edges = n_lons*(n_lats + 1) + (n_lons + 1)*n_lats
    edge_lons = np.empty((n_edges,), np.float64)
    edge_lats = np.empty((n_edges,), np.float64)
    enc_array = np.empty((n_edges, 2), np.int)
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

    print(f'lons={lons}')
    print(f'lats={lats}')
    print(f'fnc={fnc}')
    print(f'enc={enc}')

    mesh = Mesh(2, ((lons, "x"), (lats, "y")), [fnc, enc])

    mesh.add_coords(
        face_x=face_lon_coord,
        face_y=face_lat_coord,
        edge_x=edge_lon_coord,
        edge_y=edge_lat_coord,
    )
    mesh.long_name = "example mesh"
    return mesh



def _gridlike_mesh(n_lons, n_lats):
    """
    Generate a global mesh with geometry similar to a rectilinear grid.
    The resulting mesh will have n_lons cells spanning its longitudes and
    n_lats cells spanning its latitudes for a total of (n_lons * n_lats) cells.
    Note that the cells neighbouring the poles will actually be triangular while
    the rest of the cells will be rectangular.

    This function is copied from iris-esmf-regrid.
    """
    # Arrange the indices of the non-pole nodes in an array representative of their
    # latitude/longitude.
    fnc_template = np.arange((n_lats - 1) * n_lons).reshape(n_lats - 1, n_lons) + 1
    fnc_array = np.empty([n_lats, n_lons, 4])
    # Assign points in an anticlockwise orientation. From the 0 node to 1
    # longitude increases, then from 1 to 2 latitude increases, from 2 to 3
    # longitude decreases and from 3 to 0 latitude decreases.
    fnc_array[1:, :, 0] = fnc_template
    fnc_array[1:, :, 1] = np.roll(fnc_template, -1, 1)
    fnc_array[:-1, :, 2] = np.roll(fnc_template, -1, 1)
    fnc_array[:-1, :, 3] = fnc_template
    # Define the poles as single points. Note that all the cells adjacent to these
    # nodes will be topologically triangular with the pole node repeated. One of
    # these repeated pole node references will eventually be masked.
    fnc_array[0, :, :2] = 0
    num_nodes = fnc_template.max()
    fnc_array[-1, :, 2:] = num_nodes + 1
    # By convention, node references to be masked should be last in the list.
    # Since one of the pole node references will end up masked, this should be
    # moved to the end of the list of nodes.
    fnc_array[0, :, :] = np.roll(fnc_array[0, :, :], -1, -1)

    # One of the two references to the pole node are defined to be masked.
    fnc_mask = np.zeros_like(fnc_array)
    fnc_mask[0, :, -1] = 1
    fnc_mask[-1, :, -1] = 1
    fnc_ma = ma.array(fnc_array, mask=fnc_mask, dtype=int)

    # The face node connectivity is flattened to the correct dimensionality.
    fnc_ma = fnc_ma.reshape([-1, 4])

    # Describe the edge node connectivity.
    # There are n_lats * n_lons vertical edges and (n_lats - 1) * n_lons horizontal
    # edges which we arrange into an array for convenience of calculation.
    enc_array = np.empty([(n_lats * 2) - 1, n_lons, 2], dtype=int)
    # The vertical edges make up enc_array[:n_lats].
    enc_array[1:n_lats, :, 0] = fnc_template
    enc_array[: n_lats - 1, :, 1] = fnc_template
    enc_array[0, :, 0] = 0
    # The horizontal edges make up enc_array[n_lats:].
    enc_array[n_lats - 1, :, 1] = num_nodes + 1
    enc_array[n_lats:, :, 0] = fnc_template
    enc_array[n_lats:, :, 1] = np.roll(fnc_template, -1, 1)
    # The array is flattened to its proper shape of (N, 2).
    enc_array = enc_array.reshape([-1, 2])

    # Latitude and longitude values are set.
    lat_values = np.linspace(-90, 90, n_lats + 1)
    lon_values = np.linspace(-180, 180, n_lons, endpoint=False)
    # Latitude values are broadcast to arrays with the same shape as the face node
    # connectivity node references in fnc_template.
    lon_array, lat_array = np.meshgrid(lon_values, lat_values[1:-1])
    node_lats = np.empty(num_nodes + 2)
    # Note that fnc_template is created by reshaping a list of  indices. These
    # indices refer to node_lats and node_lons which are generated by reshaping
    # lat_array and lon_array. Because of the way reshaping preserves order, there
    # is a correspondance between an index in a particular position fnc_template
    # and the latitude and longitude described by lat_array and lon_array in the
    # same position.
    node_lats[1:-1] = lat_array.reshape([-1])
    # Define the latitude and longitude of the poles.
    node_lats[0] = lat_values[0]
    node_lats[-1] = lat_values[-1]
    node_lons = np.empty(num_nodes + 2)
    node_lons[1:-1] = lon_array.reshape([-1])
    node_lons[0] = 0
    node_lons[-1] = 0

    # Center Latitude and Longitude values are set.
    lon_centers = np.linspace(-180, 180, (2 * n_lons) + 1)[1::2]
    lat_centers = np.linspace(-90, 90, (2 * n_lats) + 1)[1::2]
    lon_center_array, lat_center_array = np.meshgrid(lon_centers, lat_centers)
    face_lons = lon_center_array.flatten()
    face_lats = lat_center_array.flatten()

    # Translate the mesh information into iris objects.
    fnc = Connectivity(
        fnc_ma,
        cf_role="face_node_connectivity",
        start_index=0,
    )
    enc = Connectivity(
        enc_array,
        cf_role="edge_node_connectivity",
        start_index=0,
    )
    lons = AuxCoord(node_lons, standard_name="longitude")
    lats = AuxCoord(node_lats, standard_name="latitude")

    # Create the mesh
    print(f'lons={lons}')
    print(f'lats={lats}')
    print(f'fnc={fnc}')
    print(f'enc={enc}')
    for icell in range(fnc.shape[0]):
        print(f'cell {icell} has connectivity {fnc[icell,:]}')
    mesh = Mesh(2, ((lons, "x"), (lats, "y")), [fnc, enc])

    # In order to add a mesh to a cube, face locations must be added.
    face_lon_coord = AuxCoord(face_lons, standard_name="longitude")
    face_lat_coord = AuxCoord(face_lats, standard_name="latitude")

    # Add dummy edge coords.
    dummy_points = np.zeros(enc_array.shape[0])
    edge_lon_coord = AuxCoord(dummy_points, standard_name="longitude")
    edge_lat_coord = AuxCoord(dummy_points, standard_name="latitude")

    mesh.add_coords(
        face_x=face_lon_coord,
        face_y=face_lat_coord,
        edge_x=edge_lon_coord,
        edge_y=edge_lat_coord,
    )
    mesh.long_name = "example mesh"
    return mesh


def _gridlike_mesh_cube(n_lons, n_lats, location="edge"):
    mesh = _gridlike_mesh2(n_lons, n_lats)
    mesh_coord_x, mesh_coord_y = mesh.to_MeshCoords(location)
    data = np.zeros_like(mesh_coord_x.points)
    cube = Cube(data)
    cube.add_aux_coord(mesh_coord_x, 0)
    cube.add_aux_coord(mesh_coord_y, 0)
    # cube has a mesh (cube.mesh)
    return cube

def test_cube_mesh():
    cube = _gridlike_mesh_cube(2, 3)
    assert hasattr(cube, 'shape')
    assert hasattr(cube, 'data')
    assert hasattr(cube, 'mesh')


def xtest_mesh_to_mesh_basic():
    src = _gridlike_mesh_cube(4, 5)
    tgt = _gridlike_mesh_cube(6, 3)
    src_mesh = src.mesh
    tgt_mesh = tgt.mesh

    # compute the regridding weights
    rg = mint.IrisMintRegridder(src_mesh, tgt_mesh)

    # extensive field regridding
    out_cube = rg.regrid_extensive_field(src)

    # for extensive fields, data and out_data would represent the extensive fields
    # aka flux integrals over edges
    # out_data = rg(data, data_type='extensive')

    # other option, the regridder takes vector fields on edges (u, v)
    # (Denis, u=eastward, v=northward)
    # out_u, out_v = rg((u, v), data_type='n,e components', function_space='w2') 

    # maybe if there is demand
    # other option, the regridder takes the perp component of the vector field on each edge
    # out_v = rg(v, data_type='perp component') # W2
    # out_v = rg(v, data_type='parallel component') # W1


    # res = src.regrid(tgt, MINTScheme())
    # expected = _gridlike_mesh_cube(6, 3)
    # assert res == expected


def xtest_mesh_to_mesh_streamfunction():

    # Create source cubes on unstructured meshes.
    src = _gridlike_mesh_cube(20, 10)
    tgt = _gridlike_mesh_cube(30, 20)
    rg = mint.IrisMintRegridder(src.mesh, tgt.mesh)

    rg.src.get_grid().dump('src.vtk')
    rg.src.get_grid().dump('tgt.vtk')

    # Set the edge data from a stream function [= sin(theta) + cos(theta)*cos(lambda)].
    _set_extensive_field_data_from_streamfct(src)
    _set_extensive_field_data_from_streamfct(tgt)

    # Regrid.
    result = rg.regrid_extensive_field(src)

    mint.printLogMessages()

    # Check the result.
    error = np.mean(np.fabs(result.data - tgt.data))
    print(f'error = {error}')

    #assert error < 1.e-3
    



