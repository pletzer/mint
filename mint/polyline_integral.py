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

    def build(self, grid, xyz, counterclock, periodX):
        """
        Build the flux calculator.

        :param grid: instance of Grid
        :param xyz: numpy array with npoints rows and 3 columns
        :param counterclock: orientation of the edges in the cell
                             (True=counterclockwise,
                              False=positive in xi)
        :param periodX: periodicity length in x (longitudes),
                        Set to zero if non-periodic.
        """
        MINTLIB.mnt_polylineintegral_build.argtypes = [POINTER(c_void_p), c_void_p,
                                                   c_int,
                                                   DOUBLE_ARRAY_PTR,
                                                   c_int, c_double]
        cc = 0
        if counterclock:
            cc = 1
        ier = MINTLIB.mnt_polylineintegral_build(self.obj, grid.ptr, xyz.shape[0],
                                             xyz, cc, periodX)
        if ier:
            msg = "Some target lines fall outside the grid. (Ok if these are partially outside.)"
            warning_handler(FILE, 'build', ier, detailedmsg=msg)

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
