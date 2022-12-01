from ctypes import (c_void_p, c_double, c_int, byref, POINTER,
                   c_char_p, c_size_t)
from . import MINTLIB, UNIQUE_EDGE_DATA, CELL_BY_CELL_DATA, NUM_EDGES_PER_QUAD
from . import error_handler, warning_handler
import numpy


FILE = 'regrid_edges.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class RegridEdges(object):
    """
    A class to interpolate a 1 form from a grid to another grid.
    """

    def __init__(self):
        """
        Constructor.
        """

        self.obj = byref(c_void_p())

        ier = MINTLIB.mnt_regridedges_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        ier = MINTLIB.mnt_regridedges_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setSrcGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the source grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length
                                     should be added/subtracted
                                     to make each cell as compact as
                                     possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles
                                 should be the average of the cells'
                                 longitudes

        note:: a lon-lat grid requires 0, 0 and a cubed sphere grid
               requires 1, 1
        """
        MINTLIB.mnt_regridedges_setSrcGridFlags.argtypes = [POINTER(c_void_p),
                                                        c_int, c_int]
        ier = MINTLIB.mnt_regridedges_setSrcGridFlags(self.obj,
                                                  fixLonAcrossDateline,
                                                  averageLonAtPole)
        if ier:
            error_handler(FILE, 'setSrcGridFlags', ier)

    def setDstGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the destination grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length
                                     should be added/subtracted
                                     to make each cell as compact as
                                     possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles
                                 should be the average of the cells'
                                 longitudes

        note:: a lon-lat grid requires 0, 0 and a cubed sphere grid
               requires 1, 1
        """
        MINTLIB.mnt_regridedges_setDstGridFlags.argtypes = [POINTER(c_void_p),
                                                        c_int, c_int]
        ier = MINTLIB.mnt_regridedges_setDstGridFlags(self.obj,
                                                  fixLonAcrossDateline,
                                                  averageLonAtPole)
        if ier:
            error_handler(FILE, 'setDstGridFlags', ier)

    def setSrcGrid(self, grid):
        """
        Set the source grid

        :param grid: grid
        """
        MINTLIB.mnt_regridedges_setSrcGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_regridedges_setSrcGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setSrcGrid', ier)


    def setDstGrid(self, grid):
        """
        Set the destination grid

        :param grid: grid
        """
        MINTLIB.mnt_regridedges_setDstGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_regridedges_setDstGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setDstGrid', ier)


    def loadSrcGrid(self, filename):
        """
        Load a source grid from a 2D UGRID file.

        :param filename: string in the format filename$meshname, 
                         filename should be a NetCDF file storing data in the UGRID format
                         and meshname is the name of the grid inside the file.
        """
        MINTLIB.mnt_regridedges_loadSrcGrid.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = MINTLIB.mnt_regridedges_loadSrcGrid(self.obj, fn, len(fn))

    def loadDstGrid(self, filename):
        """
        Load a destination grid from a 2D UGRID file.

        :param filename: string in the format filename$meshname,
                         filename should be a NetCDF file storing data in the UGRID format
                         and meshname is the name of the grid inside the file.
        """
        MINTLIB.mnt_regridedges_loadDstGrid.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = MINTLIB.mnt_regridedges_loadDstGrid(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'loadDstGrid', ier)

    def getNumSrcEdges(self):
        """
        Get the number of unique edges in the source grid.

        :returns number
        """
        MINTLIB.mnt_regridedges_getNumSrcEdges.argtypes = [POINTER(c_void_p)]
        MINTLIB.mnt_regridedges_getNumSrcEdges.restype = c_size_t
        n = c_size_t()
        ier = MINTLIB.mnt_regridedges_getNumSrcEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumSrcEdges', ier)
        return n.value

    def getNumDstEdges(self):
        """
        Get the number of unique edges in the destination grid.

        :returns number
        """
        MINTLIB.mnt_regridedges_getNumDstEdges.argtypes = [POINTER(c_void_p)]
        MINTLIB.mnt_regridedges_getNumDstEdges.restype = c_size_t
        n = c_size_t()
        ier = MINTLIB.mnt_regridedges_getNumDstEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumDstEdges', ier)
        return n.value

    def buildLocator(self, numCellsPerBucket=100, periodX=360., enableFolding=0):
        """
        Build the locator.

        :param numCellsPerBucket: average number of cells per bucket,
                                  performance typically improves with a higher
                                  number of cells per bucket
        :param periodX: periodicity length (set to 0 if non-periodic)
        :param enableFolding: whether (1) or not (0) to allow for |latitude| > 90
        :warning: There are cases at very coarse resolution where the regridding
                  may fail for some edges if the number of cells per bucket is too small
        """

        enableFoldingInt = 0
        if enableFolding:
            enableFoldingInt = 1

        MINTLIB.mnt_regridedges_buildLocator.argtypes = [POINTER(c_void_p), c_int,
                                                         c_double, c_int]
        ier = MINTLIB.mnt_regridedges_buildLocator(self.obj, numCellsPerBucket,
                                                   periodX, enableFoldingInt)
        if ier:
            msg = "Failed to build locator"
            warning_handler(FILE, 'buildLocator', ier, detailedmsg=msg)

    def computeWeights(self, debug=0):
        """
        Compute the regridding weights.

        :param debug: 0=no debug info, 1=print debug info, 2=save bad
                      edges in VTK file
        """
        MINTLIB.mnt_regridedges_computeWeights.argtypes = [POINTER(c_void_p), c_int]
        ier = MINTLIB.mnt_regridedges_computeWeights(self.obj, debug)
        if ier:
            msg = "Some target lines fall outside the grid. (Ok if these are partially outside.)"
            warning_handler(FILE, 'computeWeights', ier, detailedmsg=msg)

    def dumpWeights(self, filename):
        """
        Dump the weights to a file.

        :param filename: file name
        """
        MINTLIB.mnt_regridedges_dumpWeights.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = MINTLIB.mnt_regridedges_dumpWeights(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'dumpWeights', ier)

    def loadWeights(self, filename):
        """
        Load the weights from a file.

        :param filename: file name
        """
        MINTLIB.mnt_regridedges_loadWeights.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = MINTLIB.mnt_regridedges_loadWeights(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'loadWeights', ier)

    def apply(self, srcdata, dstdata, placement):
        """
        Apply the regridding weights to an edge field

        :param srcdata: contiguous array of source field data (input), see below for expected size
        :param dstdata: contiguous array of destination field data (output), see below for expected size
        :param placement: mint.CELL_BY_CELL_DATA if the data are cell by cell
                          (size of array is numCells * mint.NUM_EDGES_PER_QUAD),
                          assume unique edge Id data otherwise (size is numEdges)
        """

        MINTLIB.mnt_regridedges_apply.argtypes = [POINTER(c_void_p),
                                                  DOUBLE_ARRAY_PTR,
                                                  DOUBLE_ARRAY_PTR, c_int]
        ier = MINTLIB.mnt_regridedges_apply(self.obj, srcdata, dstdata, placement)
        if ier:
            error_handler(FILE, 'apply', ier)

    def vectorApply(self, src_u, src_v, dst_u, dst_v, fs):
        """
        Regrid a vector field

        :param src_u: x-component of the vector field on the source grid (input), size is num edges on source grid
        :param src_v: y-component of the vector field on the source grid (input), size is num edges on source grid
        :param dst_u: x-component of the vector field on the source grid (output), size is num edges on destination grid
        :param dst_v: y-component of the vector field on the source grid (output), size is num edges on destination grid
        :param fs: function space, either FUNC_SPACE_W1 or FUNC_SPACE_W2
        """

        MINTLIB.mnt_regridedges_vectorApply.argtypes = [POINTER(c_void_p),
                                                  DOUBLE_ARRAY_PTR,
                                                  DOUBLE_ARRAY_PTR,
                                                  DOUBLE_ARRAY_PTR,
                                                  DOUBLE_ARRAY_PTR,
                                                  c_int]
        ier = MINTLIB.mnt_regridedges_vectorApply(self.obj, src_u, src_v, dst_u, dst_v, fs)
        if ier:
            error_handler(FILE, 'vectorApply', ier)

