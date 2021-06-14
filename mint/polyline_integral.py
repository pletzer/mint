from ctypes import c_void_p, c_double, c_int, byref, POINTER
from . import LIB
import numpy
import logging

FILE = 'polyline_integral.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


def error_handler(filename, methodname, ier):
    msg = f'ier={ier} after calling {methodname} in {filename}!'
    logging.error(msg)
    raise RuntimeError(msg)


def warning_handler(filename, methodname, ier, detailedMsg=''):
    msg = f'ier={ier} after calling {methodname} in {filename}!\n'
    msg += detailedMsg
    logging.warning(msg)


class PolylineIntegral(object):
    """
    A class to compute line or flux integrals
    """

    def __init__(self):
        """
        Regrid edge field constructor
        """

        self.obj = byref(c_void_p())

        LIB.mnt_polylineintegral_new.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_polylineintegral_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Regrid edge field destructor
        """
        LIB.mnt_polylineintegral_del.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_polylineintegral_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def build(self, grid, xyz, counterclock, periodX):
        """
        Build the flux calculator

        :param grid: instance of Grid
        :param xyz: numpy array with npoints rows and 3 columns
        :param counterclock: orientation of the edges in the cell
                             (True=counterclockwise,
                              False=positive in xi)
        :param periodX: periodicity length in x (longitudes),
                        Set to zero if non-periodic.
        """
        LIB.mnt_polylineintegral_build.argtypes = [POINTER(c_void_p), c_void_p,
                                                   c_int,
                                                   DOUBLE_ARRAY_PTR,
                                                   c_int, c_double]
        cc = 0
        if counterclock:
            cc = 1
        ier = LIB.mnt_polylineintegral_build(self.obj, grid.ptr, xyz.shape[0],
                                             xyz, cc, periodX)
        if ier:
            msg = "Some target lines fall outside the grid."
            warning_handler(FILE, 'build', ier, detailedMsg=msg)

    def getIntegral(self, data):
        """
        Get the flux integral over the polyline

        :param data: edge field data. This array is expected to be
                     dimensioned (ncells, 4). Each value is a scalar
                     representing the integral of the field over the
                     edge. The directions of the edge are
                     (0, 0) -> (1, 0),
                     (1, 0) -> (1, 1),
                     (0, 1) -> (1, 1) and
                     (0, 0) -> (0, 1) in parametric space
        :returns the line/flux integral
        """
        LIB.mnt_polylineintegral_getIntegral.argtypes = [POINTER(c_void_p),
                                                         DOUBLE_ARRAY_PTR]
        res = c_double()
        ier = LIB.mnt_polylineintegral_getIntegral(self.obj, data, byref(res))
        if ier:
            error_handler(FILE, 'getIntegral', ier)
        return res.value
