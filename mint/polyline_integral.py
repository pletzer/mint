from ctypes import (c_void_p, c_double, c_int, byref, POINTER)
from . import MINTLIB
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

    def computeWeights(self, xyz, counterclock):
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
        if ier:
            msg = f"Failed to locate points {xyz} (ok if some fall outside the domain)"
            warning_handler(FILE, 'computeWeights', ier, detailedmsg=msg)

    def getIntegral(self, data):
        """
        Get the flux integral over the polyline.

        :param data: edge field data. This array is expected to be
                     dimensioned (numCells, 4). Each value is a scalar
                     representing the integral of the field over the
                     edge. The directions of the edges are
                     (0, 0) -> (1, 0),
                     (1, 0) -> (1, 1),
                     (0, 1) -> (1, 1) and
                     (0, 0) -> (0, 1) in parametric space
        :returns the line/flux integral
        """
        MINTLIB.mnt_polylineintegral_getIntegral.argtypes = [POINTER(c_void_p),
                                                         DOUBLE_ARRAY_PTR]
        res = c_double()
        ier = MINTLIB.mnt_polylineintegral_getIntegral(self.obj, data, byref(res))
        if ier:
            error_handler(FILE, 'getIntegral', ier)
        return res.value
