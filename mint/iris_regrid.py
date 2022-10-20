from cgitb import enable
import copy

from iris.cube import Cube
from iris.coords import DimCoord
import numpy as np

import mint

class IrisToMINTMeshAdaptor:

    def __init__(self, iris_mesh, **kwargs):
        """
        Create a MINT grid from an Iris mesh
        """

        x = iris_mesh.node_coords.node_x.points
        y = iris_mesh.node_coords.node_y.points
        num_points = x.shape[0]
        # needs to be 3d
        self.points = np.zeros((num_points, 3), np.float64)
        self.points[:, 0] = x
        self.points[:, y] = y

        self.face2node = np.array(iris_mesh.face_node_connectivity.indices_by_location(), 
            np.int64)
        self.edge2nodes = np.array(iris_mesh.edge_node_connectivity.indices_by_location(), 
            np.int64)

        # connectivity can be 1 or 0 based
        self.face2node -= iris_mesh.face_node_connectivity.start_index
        self.edge2node -= iris_mesh.edge_node_connectivity.start_index
        
        self.grid = mint.Grid()
        fixLonAcrossDateline = kwargs.get('fixLonAcrossDateline', 0)
        averageLonAtPole = kwargs.get('averageLonAtPole', 0)
        degrees = False
        if iris_mesh.node_coords.node_x.points.units == 'degrees':
            degrees = True
        self.grid.setFlagsflags(fixLonAcrossDateline, averageLonAtPole, degrees)
        self.grid.loadFromUgrid2DData(self.points, self.face2node, self.edge2node)

        self.num_faces = self.face2node.shape[0]
        self.num_edges = self.edge2nodes.shape[0]
        self.num_points = num_points

    def get_grid(self):
        return self.grid


class IrisToMINTDataAdaptor:

    def __init__(self, grid, function_space, **kwargs):
        
        self.efc = mint.ExtensiveFieldConverter()
        self.efc.setGrid(grid)

        self.getData = self.efc.getFaceData
        if function_space.lower == 'w1':
            self.getData = self.efc.getEdgeData

    def get_extensive_data(self, u, v):
        return self.getData(u, v, placement=mint.UNIQUE_EDGE_DATA)



class _DummyMintRegridder:
    def __init__(self, src_coords, tgt_coords, **kwargs):
        self.shape = tgt_coords[0].shape

    def regrid(self, data, dims, **kwargs):
        new_shape = list(data.shape)
        for dim, size in zip(dims, self.shape):
            new_shape[dim] = size
        return np.zeros(self.shape)

class _MintRegridder:

    def __init__(self, src_mesh, tgt_mesh, **kwargs):
        """
        Constructor.
        :param src_mesh: source mesh with coordinates and connectivity
        :param tgt_mesh: target mesh with coordinates and connectivity
        """

        self.src = IrisToMINTMeshAdaptor(src_mesh)
        self.tgt = IrisToMINTMeshAdaptor(tgt_mesh)

        # build the regridder
        self.regridder = mint.RegridEdges()
        self.regridder.setSrcGrid(self.src.get_grid())
        self.regridder.setDstGrid(self.tgt.get_grid())
        numCellsPerBucket = kwargs.get('numCellsPerBucket', 128)
        periodX = kwargs.get('periodX', 360.)
        enableFolding = kwargs.get('enableFolding', 0)
        self.regridder.buildLocator(numCellsPerBucket, periodX, enableFolding)

        # compute the regridding weights
        debug = kwargs.get('debug', 0)
        self.regridder.computeWeights(debug)

    def regrid(self, uv_data, function_space, dims, **kwargs):

        efc = IrisToMINTDataAdaptor(self.src.get_grid(), function_space, **kwargs)

        src_dims = uv_data[0].shape

        # last dimension is number of edges, all other dimensions are shared between
        # the src and tgt data
        dims = src_dims[:-1]

        src_num_edges = self.src.get_grid().get_num_edges()
        tgt_num_edges = self.tgt.get_grid().get_num_edges()
        tgt_num_faces = self.tgt.get_grid().get_num_faces()

        # allocate the output containers
        out_u = np.empty(dims + (tgt_num_edges), np.float64)
        out_v = np.empty(dims + (tgt_num_edges), np.float64)
        tgt_data = np.empty(dims + (tgt_num_faces*mint.NUM_EDGES_PER_QUAD,), np.float64)

        u, v = uv_data

        # assume last dimension to be edge Id
        # iterate over all the axes other than the edges
        mai = mint.MultiArrayIter(dims)
        mai.begin()
        for _ in range(mai.getNumIters()):

            inds = tuple(mai.getIndices())

            # index set for the src data on unique edges
            slab = inds + (slice(0, src_num_edges,))

            # index set 

            # TO DO transform components here

            # u and v on the UGRID
            src_data = \
                efc.getData(u.data[slab], v.data[slab], placement=mint.UNIQUE_EDGE_DATA)
            
            self.regridder.apply(src_data, tgt_data, placement=mint.CELL_BY_CELL_DATA)

            # divide by the edges' length

            # 


            mai.next()

        return (out_u, out_v)


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _make_mint_regridder(src_coords, tgt_coords, **kwargs):
    return _DummyMintRegridder(src_coords, tgt_coords, **kwargs)


def _regrid(data, dims, mint_regridder):

    # TODO: replace dummy operation with actual regridding.
    #  Actual operation should make use of mint_regridder.
    new_data = mint_regridder.regrid(data, dims)

    return new_data


def _create_cube(data, src_cube, src_dims, tgt_coords, num_tgt_dims):
    # This function is copied from iris-esmf-regrid.

    new_cube = Cube(data)

    if len(src_dims) == 2:
        grid_dim_x, grid_dim_y = src_dims
    elif len(src_dims) == 1:
        grid_dim_y = src_dims[0]
        grid_dim_x = grid_dim_y + 1
    else:
        raise ValueError(
            f"Source grid must be described by 1 or 2 dimensions, got {len(src_dims)}"
        )
    if num_tgt_dims == 1:
        grid_dim_x = grid_dim_y = min(src_dims)
    for tgt_coord, dim in zip(tgt_coords, (grid_dim_x, grid_dim_y)):
        if len(tgt_coord.shape) == 1:
            if isinstance(tgt_coord, DimCoord):
                new_cube.add_dim_coord(tgt_coord, dim)
            else:
                new_cube.add_aux_coord(tgt_coord, dim)
        else:
            new_cube.add_aux_coord(tgt_coord, (grid_dim_y, grid_dim_x))

    new_cube.metadata = copy.deepcopy(src_cube.metadata)

    def copy_coords(src_coords, add_method):
        for coord in src_coords:
            dims = src_cube.coord_dims(coord)
            if set(src_dims).intersection(set(dims)):
                continue
            offset = num_tgt_dims - len(src_dims)
            dims = [dim if dim < max(src_dims) else dim + offset for dim in dims]
            result_coord = coord.copy()
            # Add result_coord to the owner of add_method.
            add_method(result_coord, dims)

    copy_coords(src_cube.dim_coords, new_cube.add_dim_coord)
    copy_coords(src_cube.aux_coords, new_cube.add_aux_coord)

    return new_cube


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
    # TODO: properly handle the dimensionality of the target
    num_tgt_dims = 1
    return _create_cube(data, cube, dims, tgt_coords, num_tgt_dims)


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
