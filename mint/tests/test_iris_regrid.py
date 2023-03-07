import numpy as np
from pathlib import Path
import mint
from iris_utils import _u_v_cubes_from_ugrid_file, \
                       _set_vector_field_from_streamfct, \
                       _set_vector_field_from_potentialfct, \
                       _set_extensive_field_from_streamfct, \
                       _gridlike_mesh_cube

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
DEG2RAD = np.pi / 180.

def test_extfield():
    u, v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs8_wind.nc')
    efa = mint.ExtensiveFieldAdaptor()
    _set_vector_field_from_streamfct(u, v)
    flags = (1, 1, 1) # cubed-sphere
    adaptor = mint.IrisToMintMeshAdaptor(u.mesh, flags=flags)
    grid = adaptor.get_grid()
    num_cells = grid.getNumberOfCells()
    num_edges = grid.getNumberOfEdges()
    efa.setGrid(grid)
    ef = np.empty((num_edges,), np.float64)
    efa.fromVectorField(u.data, v.data, data=ef, placement=mint.UNIQUE_EDGE_DATA, fs=mint.FUNC_SPACE_W2)

def test_cs2lonlat_zt():

    tgt_u = _gridlike_mesh_cube(100, 50, height=2, time=2)
    tgt_v = _gridlike_mesh_cube(100, 50, height=2, time=2)

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind_zt.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)
    # mint.saveMeshVTK(src_u.mesh, 'cs128_w2_mesh.vtk')
    # mint.saveVectorFieldVTK(src_u, src_v, 'cs128_w2_vectors.vtk')

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)
    # grid options for a lon-lat grid
    tgt_flags = (0, 0, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags,
                                debug=3)

    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    mint.saveMeshVTK(result_u.mesh, 'lonlat_mesh.vtk')
    time_idx = 1
    height_idx = 0
    mint.saveVectorFieldVTK(result_u, result_v, f'lonlat_w2_vectors_t{time_idx}_z{height_idx}.vtk',
        extra_inds=(time_idx, height_idx))

    # check
    error = 0.5*( np.mean( np.fabs(tgt_u.data - result_u.data) ) + np.mean( np.fabs(tgt_v.data - result_v.data) ) )
    assert error < 0.01

def test_cs2lonlat():

    tgt_u = _gridlike_mesh_cube(100, 50)
    tgt_v = _gridlike_mesh_cube(100, 50)

    src_u, src_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)
    mint.saveMeshVTK(src_u.mesh, 'cs128_w2_mesh.vtk')
    mint.saveVectorFieldVTK(src_u, src_v, 'cs128_w2_vectors.vtk')

    # grid options for a cubed-sphere grid
    src_flags = (1, 1, 1)
    # grid options for a lon-lat grid
    tgt_flags = (0, 0, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags,
                                debug=3)

    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    mint.saveMeshVTK(result_u.mesh, 'lonlat_mesh.vtk')
    mint.saveVectorFieldVTK(result_u, result_v, 'lonlat_w2_vectors.vtk')

    # check
    error = 0.5*( np.mean( np.fabs(tgt_u.data - result_u.data) ) + np.mean( np.fabs(tgt_v.data - result_v.data) ) )
    assert error < 0.01

def test_lonlat2cs():

    src_u = _gridlike_mesh_cube(100, 50)
    src_v = _gridlike_mesh_cube(100, 50)

    tgt_u, tgt_v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs8_wind.nc')

    # w2
    _set_vector_field_from_streamfct(src_u, src_v)
    mint.saveMeshVTK(src_u.mesh, 'lonlat2cs_src_mesh.vtk')
    mint.saveVectorFieldVTK(src_u, src_v, 'lonlat2cs_src_w2_vectors.vtk')

    # grid options for a lonlat grid
    src_flags = (0, 0, 1)
    # grid options for a subed-sphere grid
    tgt_flags = (1, 1, 1)
    rg = mint.IrisMintRegridder(src_u.mesh, tgt_u.mesh, \
                                src_flags=src_flags, tgt_flags=tgt_flags,
                                debug=3)

    _set_vector_field_from_streamfct(tgt_u, tgt_v)

    result_u, result_v = rg.regrid_vector_cubes(src_u, src_v, fs=mint.FUNC_SPACE_W2)

    mint.saveMeshVTK(result_u.mesh, 'lonlat2cs_tgt_mesh.vtk')
    mint.saveVectorFieldVTK(result_u, result_v, 'lonlat2cs_tgt_w2_vectors.vtk')

    # check
    error = 0.5*( np.mean( np.fabs(tgt_u.data - result_u.data) ) + np.mean( np.fabs(tgt_v.data - result_v.data) ) )
    assert error < 0.01

def test_read_ugrid_file_w1():
    u_cube, v_cube = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs8_wind.nc')
    _set_vector_field_from_potentialfct(u_cube, v_cube)
    mint.saveMeshVTK(u_cube.mesh, 'cs8_w1_mesh.vtk')
    mint.saveVectorFieldVTK(u_cube, v_cube, 'cs8_w1_vectors.vtk')

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
    assert error < 0.06

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
    assert error < 0.02

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

    # Check.
    error = 0.5*(np.fabs(result_u.data - tgt_u.data).mean() + \
                 np.fabs(result_v.data - tgt_v.data).mean())
    print(f'test_cubedsphere8_to_cubedsphere8_w2 = {error}')
    assert error < 0.01

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

    assert error < 0.07

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
    assert error < 0.005

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
                # the 2 nodes are at the pole, u, v are ill-defined there so skip
                continue
            uval, vval = result_u.data[edge_id], result_v.data[edge_id]
            uxct, vxct = tgt_u.data[edge_id], tgt_v.data[edge_id]
            error += abs(uval - uxct)
            error += abs(vval - vxct)
            print(f'cell {icell} edge {ie} points {x0, y0} -> {x1, y1} u,v = {uval, vval} (exact {uxct, vxct})')
    error /= (2 * tgt_num_edges)
    print(f'error = {error}')
    assert error < 0.04
