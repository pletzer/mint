from iris.cube import Cube
import numpy as np
from mint import RegridEdges, Grid


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _make_mint_regridder(src_coords, tgt_coords, **kwargs):

    # coordinates must be edge centred
    src_x, src_y = src_coords
    assert(src_x.location == src_y.location == 'edge')
    tgt_x, tgt_y = tgt_coords
    assert(tgt_x.location == tgt_y.location == 'edge')

    # build the point arrays
    src_points = []
    src_face2nodes = []
    src_edge2nodes = []
    if len(src_x.shape) == len(src_y.shape) == 1:
        # axes
        xx, yy = np.meshgrid(src_x, src_y)
        n = np.proc(xx.shape)
        src_points = np.zeros((n, 3), np.float64)
        src_points[:, 0] =  xx[:]
        src_points[:, 1] =  yy[:]
        
    elif len(src_x.shape) == len(src_y.shape) == 2:
        # curvilinear coordinates
        n = np.proc(src_x.shape)
        src_points = np.zeros((n, 3), np.float64)
        src_points[:, 0] =  src_x[:]
        src_points[:, 1] =  src_y[:]
    else:
        # TODO mesh
        pass





    # build the src grid
    src_grid = Grid()
    fixLonAcrossDateline, averageLonAtPole, degrees = kwargs.get('src_grid_flags', (0, 0, 1))
    src_grid.setSrcGridFlags(fixLonAcrossDateline, averageLonAtPole, degrees)
    src_grid.loadFromUgrid2DData(src_xyz, src_face2nodes, src_edge2nodes)

    # build the tgt grid
    tgt_grid = Grid()
    fixLonAcrossDateline, averageLonAtPole, degrees = kwargs.get('tgt_grid_flags', (0, 0, 1))
    tgt_grid.setSrcGridFlags(fixLonAcrossDateline, averageLonAtPole, degrees)
    tgt_grid.loadFromUgrid2DData(tgt_xyz, tgt_face2nodes, tgt_edge2nodes)

    # build the regridder
    regridder = RegridEdges()
    regridder.setSrcGrid(src_grid)

    # compute the weights
    numCellsPerBucket, periodX, enableFolding = kwargs.get('numCellsPerBucket', 128), 
                                                kwargs.get('periodX', 360.0), 
                                                kwargs.get('enableFolding', 0)
    regridder.buildLocator(numCellsPerBucket, periodX, enableFolding)
    regridder.computeWeights()

    obj = dict(src_grid=src_grid, dst_grid=dst_grid, regridder=regridder, src_points=src_points, tgt_points=tgt_points)
    return obj


def _regrid(data, dims, regrid_info):
    tgt_coords, mint_regridder = regrid_info

    # TODO: replace dummy operation with actual regridding.
    #  Actual operation should make use of mint_regridder.
    out_shape = list[data.shape]
    out_shape[dims[0]] = tgt_coords[0].shape[0]
    new_data = np.oneslike(out_shape)

    return new_data


def _create_cube(data, src, src_dims, tgt_coords):
    # TODO: Add metadata to new cube.
    result_cube = Cube(data)
    return result_cube


def _prepare(src, tgt, **kwargs):
    src_coords = _get_coords(src)
    tgt_coords = _get_coords(tgt)
    mint_regridder = _make_mint_regridder(src_coords, tgt_coords, **kwargs)
    regrid_info = [tgt_coords, mint_regridder]
    return regrid_info


def _perform(cube, regrid_info):
    tgt_coords, mint_regridder = regrid_info
    dims = _get_dims(cube)
    data = _regrid(cube.data, dims, mint_regridder)
    return _create_cube(data, cube, dims, tgt_coords)


class _MINTRegridder:
    def __init__(self, src, tgt, **kwargs):
        self.regrid_info = _prepare(src, tgt, **kwargs)

    def __call__(self, src):
        return _perform(src, self.regrid_info)


class MINTScheme:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def regridder(self, src, tgt):
        return _MINTRegridder(src, tgt, **self.kwargs)
