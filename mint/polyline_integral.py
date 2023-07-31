from ctypes import (c_void_p, c_double, c_int, byref, POINTER)
from . import MINTLIB, NUM_EDGES_PER_QUAD, UNIQUE_EDGE_DATA, CELL_BY_CELL_DATA
from . import error_handler, warning_handler
import numpy


FILE = 'polyline_integral.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class PolylineIntegral(object):
    """
    A class to compute line or flux integrals.
    """

    def __init__(self):
        """
        Constructor.
        """

        self.obj = byref(c_void_p())

        MINTLIB.mnt_polylineintegral_new.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_polylineintegral_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_polylineintegral_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_polylineintegral_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setGrid(self, grid):
        """
        Set the grid.

        :param grid: instance of Grid
        """
        MINTLIB.mnt_polylineintegral_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_polylineintegral_setGrid(self.obj, grid.ptr)
        if ier:
            msg = "Failed to set the grid"
            warning_handler(FILE, 'setGrid', ier, detailedmsg=msg)

    def buildLocator(self, numCellsPerBucket=128, periodX=0.0, enableFolding=False):
        """
        Build the locator.

        :param numCellsPerBucket: average number of cells per bucket,
                                  performance typically improves with a higher
                                  number of cells per bucket
        :param periodX: periodicity length in x (longitudes),
                        set to 0 if non-periodic.
        :param enableFolding: whether (1) or not (0) to allow for |latitude| > 90
        :note: call this after setGrid
        """

        enableFoldingInt = 0
        if enableFolding:
            enableFoldingInt = 1

        MINTLIB.mnt_polylineintegral_buildLocator.argtypes = [POINTER(c_void_p),
                                                              c_int, c_double, c_int]
        ier = MINTLIB.mnt_polylineintegral_buildLocator(self.obj,
                                                        numCellsPerBucket, periodX, enableFoldingInt)
        if ier:
            msg = "Failed to build locator"
            warning_handler(FILE, 'buildLocator', ier, detailedmsg=msg)

    def computeWeights(self, xyz, counterclock=False):
        """
        Build the flux calculator.

        :param xyz: numpy array with npoints rows and 3 columns
        :param counterclock: orientation of the edges in the cell
                             (True=counterclockwise,
                              False=positive in xi)
        """
        MINTLIB.mnt_polylineintegral_computeWeights.argtypes = [POINTER(c_void_p),
                                                                c_int, DOUBLE_ARRAY_PTR,
                                                                c_int]
        cc = 0
        if counterclock:
            cc = 1
        ier = MINTLIB.mnt_polylineintegral_computeWeights(self.obj, xyz.shape[0],
                                                          xyz, cc)
        if ier == -1:
            msg = f"Need at least two points, got {xyz.shape[0]} point(s)"
            error_handler(FILE, 'computeWeights', ier, detailedmsg=msg)
        elif ier == -2:
            msg = f"Need to call buildLocator before invoking computeWeights"
            error_handler(FILE, 'computeWeights', ier, detailedmsg=msg)

    def getIntegral(self, data, placement):
        """
        Get the flux integral over the polyline.

        :param data: edge integrated field data. This array is expected 
                     to be dimensioned either (numCells, mint.NUM_EDGES_PER_QUAD)
                     if placement is mint.CELL_BY_CELL_DATA, or dimensioned
                     (numEdges,) if placement is mint.UNIQUE_EDGE_DATA. In the 
                     case where placement is mint.CELL_BY_CELL_DATA, the direction
                     of the edges are
                     (0, 0) -> (1, 0),
                     (1, 0) -> (1, 1),
                     (0, 1) -> (1, 1) and
                     (0, 0) -> (0, 1) in parametric space
        :param placement: mint.CELL_BY_CELL_DATA if the data are cell by cell or
                          mint.UNIQUE_EDGE_DATA if each edge has a unique Id. 
        :returns the line/flux integral
        """
        MINTLIB.mnt_polylineintegral_getIntegral.argtypes = [POINTER(c_void_p),
                                                             DOUBLE_ARRAY_PTR, c_int,
                                                             POINTER(c_double)]
        res = c_double()
        ier = MINTLIB.mnt_polylineintegral_getIntegral(self.obj, data, placement, byref(res))
        if ier:
            error_handler(FILE, 'getIntegral', ier)
        return res.value


    def vectorGetIntegral(self, u, v, fs):
        """
        Get the flux integral over the polyline, given an edge centred vector field.

        :param u: eastward component of the vector field, size num edges
        :param v: northward component of the vector field, size num edges
        :param fs: function space, either FUNC_SPACE_W1 or FUNC_SPACE_W2

        :returns the line/flux integral
        """
        MINTLIB.mnt_polylineintegral_vectorGetIntegral.argtypes = [POINTER(c_void_p),
                                                                   DOUBLE_ARRAY_PTR, DOUBLE_ARRAY_PTR, c_int,
                                                                   POINTER(c_double)]
        res = c_double()
        ier = MINTLIB.mnt_polylineintegral_vectorGetIntegral(self.obj, u, v, fs, byref(res))
        if ier:
            error_handler(FILE, 'vectorGetIntegral', ier)
        return res.value
