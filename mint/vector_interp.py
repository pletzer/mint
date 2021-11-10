from ctypes import (c_void_p, c_int, byref, POINTER,
                    c_double, c_size_t)
from . import MINTLIB
from . import error_handler, warning_handler
import numpy


FILE = 'vectorinterp.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class VectorInterp(object):
    """
    A class to compute the vector representation of a 1-form in 2D
    """

    def __init__(self):
        """
        Constructor.
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)
        self.numTargetPoints = 0
        self.numGridEdges = 0
        self.numGridCells = 0

        MINTLIB.mnt_vectorinterp_new.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_vectorinterp_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_vectorinterp_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_vectorinterp_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setGrid(self, grid):
        """
        Set the grid.

        :param grid: instance of Grid
        """
        self.numGridEdges = grid.getNumberOfEdges()
        self.numGridCells = grid.getNumberOfCells()
        MINTLIB.mnt_vectorinterp_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_vectorinterp_setGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setGrid', ier)

    def buildLocator(self, numCellsPerBucket=10, periodX=360.):
        """
        Build the cell locator.

        :param numCellsPerBucket: approximate number of cells per bucket
        :param periodX: periodicity in x (set to 0 if non-periodic)
        :note: call this after setGrid
        """
        MINTLIB.mnt_vectorinterp_buildLocator.argtypes = [POINTER(c_void_p),
                                                      c_int, c_double]
        ier = MINTLIB.mnt_vectorinterp_buildLocator(self.obj,
                                                numCellsPerBucket, periodX)
        if ier:
            error_handler(FILE, 'buildLocator', ier)

    def findPoints(self, targetPoints, tol2=1.e-12):
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
        MINTLIB.mnt_vectorinterp_findPoints.argtypes = [POINTER(c_void_p),
                                                    c_size_t,
                                                    DOUBLE_ARRAY_PTR,
                                                    c_double]
        self.numTargetPoints = targetPoints.shape[0]
        numBad = MINTLIB.mnt_vectorinterp_findPoints(self.obj,
                                                 self.numTargetPoints,
                                                 targetPoints, tol2)
        return numBad

    def getEdgeVectors(self, data, placement):
        """
        Get the edge vectors at given target points.

        :param data: edge data array of size number of unique edges
        :param placement: 0 if data are cell by cell (size num cells * 4), 
                          assume unique edge Id data otherwise (size num 
                          edges)
        :returns vector array of size numTargetPoints times 3
        :note: call this after invoking findPoints.
        """

        # check data size
        n = data.shape[-1]
        if placement == 0 and n != self.numGridCells * 4:
            msg = f"data has wrong size (= {n}), num cells*4 = {self.numGridCells*4}"
            ier = 10
            error_handler(FILE, 'getEdgeVectors', ier, detailedmsg=msg)
            return
        elif placement != 0 and n != self.numGridEdges:
            msg = f"data has wrong size (= {n}), num edges = {self.numGridEdges}"
            ier = 11
            error_handler(FILE, 'getEdgeVectors', ier, detailedmsg=msg)
            return

        MINTLIB.mnt_vectorinterp_getEdgeVectors.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int]
        res = numpy.zeros((self.numTargetPoints, 3), numpy.float64)
        ier = MINTLIB.mnt_vectorinterp_getEdgeVectors(self.obj, data, res,
                                                      placement)
        if ier:
            msg = "Some target lines fall outside the grid."
            warning_handler(FILE, 'getEdgeVectors', ier,
                            detailedmsg=msg)
        return res

    def getFaceVectors(self, data, placement):
        """
        Get the lateral face vectors at given target points.

        :param data: edge data array of size number of unique edges
        :param placement: 0 if data are cell by cell (size num cells * 4), 
                          assume unique edge Id data otherwise (size num 
                          edges)
        :returns vector array of size numTargetPoints times 3
        :note: call this after invoking findPoints.
        """

        # check data size
        n = data.shape[-1]
        if placement == 0 and n != self.numGridCells * 4:
            msg = f"data has wrong size (= {n}), num cells*4 = {self.numGridCells*4}"
            ier = 10
            error_handler(FILE, 'getEdgeVectors', ier, detailedmsg=msg)
            return
        elif placement != 0 and n != self.numGridEdges:
            msg = f"data has wrong size (= {n}), num edges = {self.numGridEdges}"
            ier = 11
            error_handler(FILE, 'getEdgeVectors', ier, detailedmsg=msg)
            return

        MINTLIB.mnt_vectorinterp_getFaceVectors.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int]
        res = numpy.zeros((self.numTargetPoints, 3), numpy.float64)
        ier = MINTLIB.mnt_vectorinterp_getFaceVectors(self.obj, data, res,
                                                      placement)
        if ier:
            msg = "Some target lines fall outside the grid."
            warning_handler(FILE, 'getFaceVectors', ier,
                            detailedmsg=msg)
        return res
