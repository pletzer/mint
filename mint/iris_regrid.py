from iris.cube import Cube
from iris.coords import AuxCoord
import numpy as np

import mint


class IrisMintRegridder:

    def __init__(self, src_mesh, tgt_mesh, src_flags, tgt_flags, **kwargs):
        """
        Constructor.
        :param src_mesh: source iris mesh with coordinates and connectivity
        :param tgt_mesh: target iris mesh with coordinates and connectivity
        :param src_flags: flags to pass to the source grid. Example (0, 0, 1) for a regular grid 
                      and (1, 1, 1) for a cubed-sphere grid. 
        :param tgt_flags: flags to pass to the target grid. Example (0, 0, 1) for a regular grid 
                      and (1, 1, 1) for a cubed-sphere grid. 
        """

        self.tgt_mesh = tgt_mesh

        self.src = mint.IrisToMintMeshAdaptor(src_mesh, flags=src_flags)
        self.tgt = mint.IrisToMintMeshAdaptor(tgt_mesh, flags=tgt_flags)

        self.src_num_edges = self.src.get_grid().getNumberOfEdges()
        self.tgt_num_edges = self.tgt.get_grid().getNumberOfEdges()

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


    def regrid_vector_cubes(self, u_cube, v_cube, fs, **kwargs):
        """
        Regrid a vector field
        :param u_cube: eastward component of the vector fields defined on horizontal edges
        :param v_cube: northward component of the vector fields defined on horizontal edges
        :param fs: function space, either 1 (for W1/edge) or 2 (for W2/face)
        :returns (u, v) cubes
        """

        if not isinstance(u_cube.coords()[-1], AuxCoord) or not isinstance(v_cube.coords()[-1], AuxCoord):
            msg = f'Last coordinate must be of type AuxCoord'
            raise ValueError(msg)

        if np.any(u_cube.shape != v_cube.shape):
            msg = f'Source u,v cubes must have the same dimensions'
            raise ValueError(msg)


        # Dimensions other than horizontal
        dims = u_cube.shape[:-1] # last dimension is assumed to be the number of edges

        # Allocate.
        tgt_u_data = np.empty(dims + (self.tgt_num_edges,), np.float64)
        tgt_v_data = np.empty(dims + (self.tgt_num_edges,), np.float64)
        
        # Iterate over the dimensions other than horizontal.
        mai = mint.MultiArrayIter(dims)
        mai.begin()
        for _ in range(mai.getNumIters()):

            inds = tuple(mai.getIndices())

            src_slab = inds + (slice(0, self.src_num_edges),)
            tgt_slab = inds + (slice(0, self.tgt_num_edges),)

            # Regrid the vector fields.
            self.regridder.vectorApply(u_cube.data[src_slab], v_cube.data[src_slab], \
                                       tgt_u_data[tgt_slab], tgt_v_data[tgt_slab], fs)
                        
            mai.next()

        # Build the cubes.
        tgt_mesh_coord_x, tgt_mesh_coord_y = self.tgt_mesh.to_MeshCoords("edge")

        out_u_cube = Cube(tgt_u_data)
        out_v_cube = Cube(tgt_v_data)

        i = 0
        for coord in u_cube.coords(dim_coords=True):
            out_u_cube.add_dim_coord(coord, i)
            out_v_cube.add_dim_coord(coord, i)
            i += 1

        n = len(u_cube.shape)
        nm1 = n - 1
        out_u_cube.add_aux_coord(tgt_mesh_coord_x, nm1)
        out_u_cube.add_aux_coord(tgt_mesh_coord_y, nm1)

        return (out_u_cube, out_v_cube)


    def regrid_extensive_cube(self, cube, **kwargs):
        """
        Regrid the extensive cube data
        :param cube: source cube on Mesh, of size of data is ..., num edges
        :returns a new cube on the target mesh
        """

        tgt_data = self.regrid_extensive_data(cube.data, **kwargs)
    
        # build the cube
        out_cube = Cube(tgt_data)
        tgt_mesh_coord_x, tgt_mesh_coord_y = self.tgt_mesh.to_MeshCoords("edge")
        out_cube.add_aux_coord(tgt_mesh_coord_x, 0)
        out_cube.add_aux_coord(tgt_mesh_coord_y, 0)

        return out_cube


    def regrid_extensive_data(self, data, **kwargs):
        """
        Regrid the extensive data
        :param data: source data of size ..., num edges
        :returns regridded data on the target mesh
        """

        # all the dimensions other than horizontal. Assuming 
        # that the last dimension is the number of edges
        dims = data.shape[:-1]
        
        tgt_data = np.empty(dims + (self.tgt_num_edges,), np.float64)

        mai = mint.MultiArrayIter(dims)
        mai.begin()
        for _ in range(mai.getNumIters()):

            inds = tuple(mai.getIndices())

            src_slab = inds + (slice(0, self.src_num_edges),)
            tgt_slab = inds + (slice(0, self.src_num_edges),)

            src_d = data[src_slab]
            tgt_d = tgt_data[tgt_slab]

            self.regridder.apply(src_d, tgt_d, placement=mint.UNIQUE_EDGE_DATA)

            mai.next()

        return tgt_data


def _get_dims(cube):
    # This function may be made more sophisticated to accomodate
    # data not located on a mesh.
    dims = [cube.mesh_dim()]
    return dims


def _get_coords(cube):
    coords = (cube.coord(axis="X"), cube.coord(axis="Y"))
    return coords


def _make_mint_regridder(src_mesh, tgt_mesh, **kwargs):
    return _MintRegridder(src_mesh, tgt_mesh, **kwargs)


def _regrid(data, dims, mint_regridder):

    new_data = mint_regridder.regrid(data, dims)
    return new_data


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
    return mint.createIrisCube(data, cube, dims, tgt_coords, num_tgt_dims)


class _MINTRegridder:
    def __init__(self, src, tgt, **kwargs):
        """
        src and tgt are cubes
        """
        self.regrid_info = _prepare(src, tgt, **kwargs)

    def __call__(self, src):
        """
        src is a cube
        returns a cube
        """
        return _perform(src, self.regrid_info)


class MINTScheme:
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def regridder(self, src, tgt):
        return _MINTRegridder(src, tgt, **self.kwargs)
