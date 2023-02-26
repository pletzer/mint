from iris.cube import Cube
from iris.coords import DimCoord, AuxCoord
import numpy as np

import mint


class IrisMintFlux:

    def __init__(self, src_mesh, src_flags, tgt_line, **kwargs):
        """
        Constructor.
        :param src_mesh: source iris mesh with coordinates and connectivity
        :param src_flags: flags to pass to the source grid. Example (0, 0, 1) for a regular grid 
                      and (1, 1, 1) for a cubed-sphere grid.
        :param kwargs: additional arguments, eg numCellsPerBucket, periodX, enableFolding...
        :param tgt_line: array of points [(x0,y0), (x1, y1), ...] where x0, y0 ... are longitude-latitude pairs
        """

        self.src = mint.IrisToMintMeshAdaptor(src_mesh, flags=src_flags)

        self.src_num_edges = self.src.get_grid().getNumberOfEdges()

        # build the flux calculator
        self.flux_calc = mint.PolylineIntegral()
        self.flux_calc.setGrid(self.src.get_grid())

        numCellsPerBucket = kwargs.get('numCellsPerBucket', 128)
        periodX = kwargs.get('periodX', 360.)
        enableFolding = kwargs.get('enableFolding', 0)
        self.flux_calc.buildLocator(numCellsPerBucket, periodX, enableFolding)

        # compute the interpolation weights
        if hasattr(tgt_line, 'shape') and len(tgt_line.shape) == 2 and tgt_line.shape[1] == 3:
            xyz = tgt_line
        else:
            xyz = np.array([(xy[0], xy[1], 0.) for xy in tgt_line])
        self.flux_calc.computeWeights(xyz)


    def evaluate_from_vector(self, u_cube, v_cube, fs):

        """
        Evaluate the flux integral from vector components.
        :param u_cube: eastward component of the vector fields defined on horizontal edges
        :param v_cube: northward component of the vector fields defined on horizontal edges
        :param fs: function space, either 1 (for W1/edge) or 2 (for W2/face)
        :returns flux array
        """

        if not isinstance(u_cube.coords()[-1], AuxCoord) or not isinstance(v_cube.coords()[-1], AuxCoord):
            msg = f'Last coordinate must be of type AuxCoord'
            raise ValueError(msg)

        if np.any(u_cube.shape != v_cube.shape):
            msg = f'Source u,v cubes must have the same dimensions'
            raise ValueError(msg)


        if not isinstance(u_cube.coords()[-1], AuxCoord) or not isinstance(v_cube.coords()[-1], AuxCoord):
            msg = f'Last coordinate must be of type AuxCoord'
            raise ValueError(msg)

        if np.any(u_cube.shape != v_cube.shape):
            msg = f'Source u,v cubes must have the same dimensions'
            raise ValueError(msg)


        # Dimensions other than horizontal
        dims = u_cube.shape[:-1] # last dimension is assumed to be the number of edges

        # Allocate.
        res = np.empty(dims, np.float64)
        
        # Iterate over the dimensions other than horizontal.
        mai = mint.MultiArrayIter(dims)
        mai.begin()
        for _ in range(mai.getNumIters()):

            inds = tuple(mai.getIndices())

            src_slab = inds + (slice(0, self.src_num_edges),)

            # Evaluate the flux.
            res[inds] = self.flux_calc.vectorGetIntegral(u_cube.data[src_slab], v_cube.data[src_slab], \
                                                        fs=fs)
                        
            mai.next()

        return res


    def evaluate_from_extensive_cube(self, cube, **kwargs):
        """
        Evaluate flux from extensive cube data.
        :param cube: source cube on Mesh, of size of data is ..., num edges
        :returns flux array
        """

        if not isinstance(cube.coords()[-1], AuxCoord):
            msg = f'Last coordinate must be of type AuxCoord'
            raise ValueError(msg)

        # Dimensions other than horizontal.
        dims = cube.shape[:-1] # last dimension is assumed to be the number of edges

        # Allocate.
        res = np.empty(dims, np.float64)
        
        # Iterate over the dimensions other than horizontal.
        mai = mint.MultiArrayIter(dims)
        mai.begin()
        for _ in range(mai.getNumIters()):

            inds = tuple(mai.getIndices())

            src_slab = inds + (slice(0, self.src_num_edges),)

            # Evaluate the flux.
            res[inds] = self.flux_calc.getIntegral(cube.data[src_slab],\
                                                  placement=mint.UNIQUE_EDGE_DATA)
                        
            mai.next()

        return res

