from ctypes import (c_void_p, c_int, byref, POINTER,
                    c_double, c_size_t)
from . import LIB
import numpy
import logging


def error_handler(filename, methodname, ier):
    msg = f'ier={ier} after calling {methodname} in {filename}!'
    logging.error(msg)
    raise RuntimeError(msg)


FILE = 'vectorinterp.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class VectorInterp(object):
    """
    A class to compute the vector representation of a 1-form in 2D
    """

    def __init__(self):
        """
        Vector interpolator constructor.
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)
        self.numTargetPoints = 0

        LIB.mnt_vectorinterp_new.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_vectorinterp_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Vector interpolator destructor.
        """
        LIB.mnt_vectorinterp_del.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_vectorinterp_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setGrid(self, grid):
        """
        Set the grid.

        :param grid: instance of a Grid
        """
        LIB.mnt_vectorinterp_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = LIB.mnt_vectorinterp_setGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setGrid', ier)

    def buildLocator(self, numCellsPerBucket, periodX):
        """
        Build the cell locator.

        :param numCellsPerBucket: approximate number of cells per bucket
        :param periodX: periodicity in x (set to 0 if non-periodic)
        :note: call this after setGrid
        """
        LIB.mnt_vectorinterp_buildLocator.argtypes = [POINTER(c_void_p),
                                                      c_int, c_double]
        ier = LIB.mnt_vectorinterp_buildLocator(self.obj,
                                                numCellsPerBucket, periodX)
        if ier:
            error_handler(FILE, 'buildLocator', ier)

    def findPoints(self, targetPoints, tol2):
        """
        Find the cells containing the target points.

        :param targetPoints: array of size numPoints times 3
        :param tol2: tolerance in the square of the distance
        :returns the number of points outside the domain
        """
        if len(targetPoints.shape) != 2:
            msg = 'ERROR: targetPoints should have dims (numTargetPoints, 3)'\
                  f', got {targetPoints.shape}'
            raise RuntimeError(msg)
        if targetPoints.shape[-1] != 3:
            msg = "ERROR: targetPoints' last dimension should be 3,"\
                   f" got {targetPoints.shape[-1]}"
            raise RuntimeError(msg)
        LIB.mnt_vectorinterp_findPoints.argtypes = [POINTER(c_void_p),
                                                    c_size_t,
                                                    DOUBLE_ARRAY_PTR,
                                                    c_double]
        self.numTargetPoints = targetPoints.shape[0]
        numBad = LIB.mnt_vectorinterp_findPoints(self.obj,
                                                 self.numTargetPoints,
                                                 targetPoints, tol2)
        return numBad

    def getVectors(self, data):
        """
        Get the vectors at the target points.

        :param data: edge data array of size numCells times 4
        :returns vector array of size numTargetPoints times 3
        :note: call this after invoking findPoints
        """
        LIB.mnt_vectorinterp_getVectors.argtypes = [POINTER(c_void_p),
                                                    DOUBLE_ARRAY_PTR,
                                                    DOUBLE_ARRAY_PTR]
        res = numpy.zeros((self.numTargetPoints, 3), numpy.float64)
        ier = LIB.mnt_vectorinterp_getVectors(self.obj, data, res)
        if ier:
            error_handler(FILE, 'getVectors', ier)
        return res
