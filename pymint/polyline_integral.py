from ctypes import c_void_p, c_double, c_int, byref, POINTER, c_char_p, c_size_t, c_longlong
from . import LIB
import numpy

def error_handler(filename, methodname, ier):
    raise RuntimeError(f'ERROR ier={ier} after calling {methodname} in {filename}!')

FILE = 'polyline_integral.py'

class PolylineIntegral(object):


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


    def setGrid(self, grid):
        """
        Set the grid

        :param grid: instance of Grid
        """
        LIB.mnt_polylineintegral_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = LIB.mnt_polylineintegral_setGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setGrid', ier)


    def setPolyline(self, xyz):
        """
        Set the target polyline

        :param xyz: numpy array with npoints rows and 3 columns
        """
        LIB.mnt_polylineintegral_setPolyline.argtypes = [POINTER(c_void_p), c_int,
                                             numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        ier = LIB.mnt_polylineintegral_setPolyline(self.obj, xyz.shape[0], xyz)
        if ier:
            error_handler(FILE, 'setPolyline', ier)


    def getIntegral(self, data):
        """
        Get the flux integral over the polyline

        :param data: edge field data. This array is expected to be dimensioned (ncells, 4). 
                     Each data value is a scalar representing the integral of the field over 
                     the edge. The directions of the edge are 
                     (0, 0) -> (1, 0), 
                     (1, 0) -> (1, 1), 
                     (0, 1) -> (1, 1) and 
                     (0, 0) -> (0, 1) in parametric space
        :returns the line/flux integral
        """
        LIB.mnt_polylineintegral_getIntegral.argtypes = [POINTER(c_void_p),
                                             numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        LIB.mnt_polylineintegral_getIntegral.restype = c_double
        res = c_double()
        ier = LIB.mnt_polylineintegral_getIntegral(self.obj, data, byref(res))
        if ier:
            error_handler(FILE, 'getIntegral', ier)
        return res.value

