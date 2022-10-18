# from iris.cube import Cube
from wsgiref import validate
import numpy as np
from mint import NUM_EDGES_PER_QUAD, NUM_VERTS_PER_QUAD, UNIQUE_EDGE_DATA, CELL_BY_CELL_DATA, RegridEdges, Grid, regrid_edges
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
    """
    Apply the regridding weights to the source field
    :param uv_data: (u, v) velocity components on the source grid, assumed to be Arakawa C
    :returns (u, v) on the target grid
    """

    # earth's radius
    planet_radius = kwargs.get('A', 6371e3)

    # compute the extensive field from the u, v vector fields
    u, v = uv_data[0].copy(), uv_data[1].copy()
    u.units, v.units = uv_data[0].units, uv_data[1].units

    #
    # assume the (u, v) components to be a velocity field in m/s
    # for the time being
    #

    deg2rad = np.pi/180.
    if (u.units == 'm s-1') and (v.units == 'm s-1'):
        # velocity, transform to rad/s
        yy = regrid_info['src_grid']['coords'][1]
        # lat at v points
        yyv = 0.5*(yy[:,:-1] + yy[:, 1:])
        # u on u points
        u /= planet_radius # Arakawa C
        # v on v points
        v /= (planet_radius * np.cos(deg2rad*yyv)) # Arakawa C
    else:
        raise RuntimeError("unsupported vector field units -- don't know to transform to deg units")

    # convert the (u, v) to be cell by cell
    src_num_cells = v.shape[1] * u.shape[0]

    src_u_cell_by_cell = np.zeros( (src_num_cells, NUM_EDGES_PER_QUAD), np.float64 )
    src_u_cell_by_cell[:, 1] = u[:, 1:].ravel()  # Arakawa C
    src_u_cell_by_cell[:, 3] = u[:, :-1].ravel() # Arakawa C

    src_v_cell_by_cell = np.zeros( (src_num_cells, NUM_EDGES_PER_QUAD), np.float64 )
    src_v_cell_by_cell[:, 0] = v[:-1, :].ravel() # Arakawa C
    src_v_cell_by_cell[:, 2] = v[1:, :].ravel()  # Arakawa C

    # flatten
    src_u_cell_by_cell = src_u_cell_by_cell.ravel()
    src_v_cell_by_cell = src_v_cell_by_cell.ravel()

    # compute the edge integrals
    efc = ExtensiveFieldConverter()
    efc.setGrid(regrid_info['src_grid']['grid'])
    src_edge_data = efc.getFaceData(src_u_cell_by_cell, src_v_cell_by_cell, placement=CELL_BY_CELL_DATA)

    tgt_num_cells = regrid_info['tgt_grid']['grid'].getNumberOfCells()

    tgt_edge_data = np.empty((tgt_num_cells, NUM_EDGES_PER_QUAD), np.float64)
    
    # apply the interpolation weights. tgt_edge_data holds the result
    regrid_info['regridder'].apply(src_edge_data, tgt_edge_data, placement=CELL_BY_CELL_DATA)

    # convert the line integrals to vector field components
    tgt_xx, tgt_yy = regrid_info['tgt_grid']['coords']
    tgt_ny1, tgt_nx1 = tgt_xx.shape
    tgt_nx = tgt_nx1 - 1
    tgt_ny = tgt_ny1 - 1
    new_u = np.empty((tgt_ny, tgt_nx1), dtype=uv_data[0].dtype)  # Arakawa C
    new_v = np.empty((tgt_ny1, tgt_nx), dtype=uv_data[1].dtype)  # Arakawa C

    tgt_edge_data = tgt_edge_data.reshape((tgt_num_cells, NUM_EDGES_PER_QUAD))

    # set u on the left side of the cells
    new_u[:, :-1] = tgt_edge_data[:, 3].reshape((tgt_ny, tgt_nx))  # Arakawa C
    # set u on the right side of the cells
    new_u[:, 1:] = tgt_edge_data[:, 1].reshape((tgt_ny, tgt_nx))  # Arakawa C

    # set v on the bottom of the cells
    new_v[:-1, :] = tgt_edge_data[:, 0].reshape((tgt_ny, tgt_nx))  # Arakawa C
    # set v on the top of the cells
    new_v[1:, :] = tgt_edge_data[:, 2].reshape((tgt_ny, tgt_nx))  # Arakawa C

    # divide by edge lengths to recover vectors from fluxes. new_u and new_v are the 
    # vector components perpendicular to the edge!
    new_u /= (tgt_yy[1:, :] - tgt_yy[:-1, :]) # Arakawa C
    new_v /= (tgt_xx[:, 1:] - tgt_xx[:, :-1]) # Arakawa C

    # convert to m/s
    new_u *= planet_radius
    new_v *= planet_radius * np.cos(0.5*(tgt_yy[:, :-1] + tgt_yy[:, 1:])*deg2rad)


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
        return _regrid(src, self.regrid_info)
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
