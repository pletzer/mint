# from iris.cube import Cube
import numpy as np
from mint import NUM_VERTS_PER_QUAD, UNIQUE_EDGE_DATA, RegridEdges, Grid, regrid_edges
from mint.extensive_field_converter import ExtensiveFieldConverter


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _build_cell_by_cell_xyz(coords):
    """
    Build the cell be cell coordinate arrays
    :param src_coords: source grid coordinates (src_x, src_y), src_x and src_y have dimensions (src_ny, src_nx)
    :param tgt_coords: target grid coordinates (tgt_x, tgt_y)
    """
    x, y = coords

    ny1, nx1 = x.shape
    ny, nx = ny1 - 1, nx1 - 1
    num_cells = ny * nx

    # allocate
    xyz = np.zeros((num_cells, NUM_VERTS_PER_QUAD, 3), np.float64)

    # vertex ordering
    #  3......2
    #  :      :
    #  0......1

    # set the lons
    xyz[:, 0, 0] = x[:-1, :-1].ravel()
    xyz[:, 1, 0] = x[:-1, 1:].ravel()
    xyz[:, 2, 0] = x[1:, 1:].ravel()
    xyz[:, 3, 0] = x[1:, :-1].ravel()

    # set the lats
    xyz[:, 0, 1] = y[:-1, :-1].ravel()
    xyz[:, 1, 1] = y[:-1, 1:].ravel()
    xyz[:, 2, 1] = y[1:, 1:].ravel()
    xyz[:, 3, 1] = y[1:, :-1].ravel()
    
    # elevs are zero

    return xyz


def _build_grid(coords, **kwargs):

    grid = Grid()

    fixLonAcrossDateline, averageLonAtPole, degrees = kwargs.get('flags', (0, 0, 1))
    grid.setFlags(fixLonAcrossDateline, averageLonAtPole, degrees)

    xyz = _build_cell_by_cell_xyz(coords)
    grid.setPoints(xyz)

    obj = dict(grid=grid, xyz=xyz, coords=coords)
    return obj


def _make_mint_regridder(src_coords, tgt_coords, **kwargs):
    """
    Make a MINT regridder
    :param src_coords: source grid coordinates (src_x, src_y), src_x and src_y have dimensions (src_ny, src_nx)
    :param tgt_coords: target grid coordinates (tgt_x, tgt_y)
    """
    # get all the flags for the src and tgt grids
    src_kwargs = kwargs.get('src_grid', {})
    tgt_kwargs = kwargs.get('tgt_grid', {})

    # build the src and tgt grid objects
    src_grid_obj = _build_grid(src_coords, **src_kwargs)
    tgt_grid_obj = _build_grid(tgt_coords, **tgt_kwargs)

    # build the regridder
    regridder = RegridEdges()
    regridder.setSrcGrid(src_grid_obj['grid'])
    regridder.setDstGrid(tgt_grid_obj['grid'])

    # compute the weights
    numCellsPerBucket = kwargs.get('numCellsPerBucket', 128)
    periodX = kwargs.get('periodX', 360.0)
    enableFolding = kwargs.get('enableFolding', 0)                     
    print(f'numCellsPerBucket={numCellsPerBucket} periodX={periodX} enableFolding={enableFolding}')
    regridder.buildLocator(numCellsPerBucket, periodX, enableFolding)
    regridder.computeWeights()

    regrid_info = dict(src_grid=src_grid_obj, tgt_grid=tgt_grid_obj, 
               regridder=regridder)
    return regrid_info


def _regrid(uv_data, regrid_info, **kwargs):

    # compute the extensive field from the u, v vector fields
    u, v = uv_data

    rad2deg = 180./np.pi
    deg2rad = np.pi/180.
    if u.units == v.units == 'm s-1':
        # velocity, transform to rad/s
        yy = regrid_info['src_grid']['coords'][1]
        # lat at u points
        yyu = 0.5*(yy[1:,:] + yy[:-1, :])
        # earth's radius
        planet_radius = kwargs.get('A', 6371e3)
        # u on u points
        u *= rad2deg / (planet_radius * np.cos(deg2rad*yyu))
        # v on v points
        v *= rad2deg / planet_radius
    else:
        raise RuntimeError("unsupported vector field units -- don't know to transform to deg units")

    src_num_u_egdes = np.prod(u.shape)
    src_num_v_edges = np.prod(v.shape)
    src_num_edges = src_num_u_egdes + src_num_v_edges

    # the first src_num_u_edges are u edges, then followed by src_num_v_edges
    src_vx_edges = np.zeros((src_num_edges,), np.float64)
    src_vx_edges[:num_u_edges] = u
    src_vy_edges = np.zeros((src_num_edges,), np.float64)
    src_vy_edges[num_u_edges:] = v

    efc = ExtensiveFieldConverter()
    efc.setGrid(regrid_info['src_grid_obj']['grid'])
    src_edge_data = efc.getFaceData(src_vx_edges, src_vy_edges, placement=mint.UNIQUE_EDGE_DATA)

    mint_regridder = regrid_info['regridder']

    tgt_xx, tgt_yy = regrid_info['tgt_grid_obj']['coords']
    tgt_num_u_edges = np.prod(tgt_xx[1:,:].shape)
    tgt_num_v_edges = np.prod(tgt_xx[:,1:].shape)
    tgt_num_edges = tgt_num_u_edges + tgt_num_v_edges
    tgt_edge_data = np.empty((tgt_num_edges,), np.float64)
    
    mint_regridder.apply(src_edge_data, tgt_edge_data, placement=UNIQUE_EDGE_DATA)

    tgt_ny1, tgt_nx1 = tgt_xx.shape
    tgt_nx = tgt_nx1 - 1
    tgt_ny = tgt_ny1 - 1
    new_u = np.empty((tgt_ny, tgt_nx1), dtype=uv_data[0].dtype())  # Arakawa C
    new_v = np.empty((tgt_ny1, tgt_nx), dtype=uv_data[1].dtype())  # Arakawa C

    # divide by edge lengths to recover vectors from fluxes
    new_u[:] = tgt_edge_data[:tgt_num_u_edges] / (tgt_yy[1:, :] - tgt_yy[:-1, :])
    new_v[:] = tgt_edge_data[tgt_num_u_edges:] / (tgt_xx[:, 1:] - tgt_xx[:, :-1])

    return (new_u, new_v)


def _create_cube(data, src, src_dims, tgt_coords):
    # TODO: Add metadata to new cube.
    result_cube = Cube(data)
    return result_cube


def _prepare(src, tgt, **kwargs):
    regrid_info = _make_mint_regridder(src, tgt, **kwargs)
    return regrid_info


# def _perform(uv, regrid_info):
#     data = _regrid(cube.data, mint_regridder)
#     return _create_cube(data, cube, dims, tgt_coords)


class _MINTRegridder:
    def __init__(self, src, tgt, **kwargs):
        """
        Create a MINT regridder
        :param src: source grid coordinates (src_x, src_y), src_x and src_y have dimensions (src_ny, src_nx)
        :param tgt: target grid coordinates (tgt_x, tgt_y)
        """
        self.regrid_info = _prepare(src, tgt, **kwargs)

    def __call__(self, src):
        """
        Apply the regridding weights
        :param src: (src_u, src_v) arrays on source grid edges
        :returns (tgt_u, tgt_v) arrays on target grid edges
        """
        return _regrid(src, self.regrid_info, **kwargs)
        # return _perform(src, self.regrid_info)


class MINTScheme:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def regridder(self, src, tgt):
        """
        Get a regridder instance
        :param src: source grid coordinates (src_x, src_y), src_x and src_y have dimensions (src_ny, src_nx)
        :param tgt: target grid coordinates (tgt_x, tgt_y)
        :returns a regridder instance
        """
        return _MINTRegridder(src, tgt, **self.kwargs)
