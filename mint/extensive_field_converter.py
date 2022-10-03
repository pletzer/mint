from ctypes import (c_void_p, c_int, byref, POINTER,
                    c_double, c_size_t)
from . import MINTLIB, UNIQUE_EDGE_DATA, CELL_BY_CELL_DATA, NUM_EDGES_PER_QUAD
from . import error_handler, warning_handler
import numpy


FILE = 'extensive_field_converter.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class ExtensiveFieldConverter(object):
    """
    A class to convert intensive fields to extensive fields
    """

    def __init__(self):
        """
        Constructor.
        """

        self.numGridEdges = 0
        self.numGridCells = 0
        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        MINTLIB.mnt_extensivefieldconverter_new.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_extensivefieldconverter_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_extensivefieldconverter_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_extensivefieldconverter_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setGrid(self, grid):
        """
        Set the grid.

        :param grid: instance of Grid
        """
        self.numGridEdges = grid.getNumberOfEdges()
        self.numGridCells = grid.getNumberOfCells()
        MINTLIB.mnt_extensivefieldconverter_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_extensivefieldconverter_setGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setGrid', ier)

    def getEdgeData(self, vx, vy, placement):
        """
        Get the edge integrated data
        :param vx: x component of vectors on edges, array (see placement argument below)
        :param vy: y component of vectors on edges, array (see placement argument below)
        :param placement: mint.CELL_BY_CELL_DATA if vx and vy are cell by cell (size num cells * mint.NUM_EDGES_PER_QUAD),
                          vx and vy are assumed to be on unique edges otherwise (size num edges)
        :returns edge integrated data, array size is numCells * mint.NUM_EDGES_PER_QUAD
        :note: vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
               velocities in m/s and the coordinates are in degrees then one needs transform 
               vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
               (A is the planet radius).
        """

        # check the data size
        n0 = self.numGridCells * NUM_EDGES_PER_QUAD
        if placement == CELL_BY_CELL_DATA and (vx.shape[0] != n0 or vy.shape[0] != n0):
            msg = f"vx, vy have wrong size (= {vx.shape[0]}, {vy.shape[0]}), num cells*NUM_EDGES_PER_QUAD = {n0}"
            ier = 10
            error_handler(FILE, 'getEdgeData', ier, detailedmsg=msg)
            return
        elif placement != CELL_BY_CELL_DATA and (vx.shape[0] != self.numGridEdges or vy.shape[0] != self.numGridEdges):
            msg = f"vx, vy have wrong size (= {vx.shape[0]}, {vy.shape[0]}), num edges = {self.numGridEdges}"
            ier = 11
            error_handler(FILE, 'getEdgeData', ier, detailedmsg=msg)
            return

        res = numpy.empty(n0, numpy.float64)
        MINTLIB.mnt_extensivefieldconverter_getEdgeData.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int,
                                                        DOUBLE_ARRAY_PTR]
        ier = MINTLIB.mnt_extensivefieldconverter_getEdgeData(self.obj, vx, vy, placement, res)
        if ier:
            msg = "An error occurred after calling MINTLIB.mnt_extensivefieldconverter_getEdgeData."
            warning_handler(FILE, 'getEdgeData', ier,
                            detailedmsg=msg)
        return res

    def getFaceData(self, vx, vy, placement):
        """
        Get the face integrated data
        :param vx: x component of vectors on edges, array (see placement argument below)
        :param vy: y component of vectors on edges, array (see placement argument below) 
        :param placement: mint.CELL_BY_CELL_DATA if vx and vy are cell by cell (size num cells * mint.NUM_EDGES_PER_QUAD),
                          vx and vy are assumed to be on unique edges otherwise (size num edges)
        :returns face integrated data, array size is numCells * mint.NUM_EDGES_PER_QUAD
        :note: vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
               velocities in m/s and the coordinates are in degrees then one needs transform 
               vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
               (A is the planet radius).
        """

        # check the data size
        n0 = self.numGridCells * NUM_EDGES_PER_QUAD
        if placement == CELL_BY_CELL_DATA and (vx.shape[0] != n0 or vy.shape[0] != n0):
            msg = f"vx, vy have wrong size (= {vx.shape[0]}, {vy.shape[0]}), num cells*NUM_EDGES_PER_QUAD = {n0}"
            ier = 10
            error_handler(FILE, 'getEdgeData', ier, detailedmsg=msg)
            return
        elif placement != CELL_BY_CELL_DATA and (vx.shape[0] != self.numGridEdges or vy.shape[0] != self.numGridEdges):
            msg = f"vx, vy have wrong size (= {vx.shape[0]}, {vy.shape[0]}), num edges = {self.numGridEdges}"
            ier = 11
            error_handler(FILE, 'getFaceData', ier, detailedmsg=msg)
            return

        res = numpy.empty(n0, numpy.float64)
        MINTLIB.mnt_extensivefieldconverter_getFaceData.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int,
                                                        DOUBLE_ARRAY_PTR]
        ier = MINTLIB.mnt_extensivefieldconverter_getFaceData(self.obj, vx, vy, placement, res)
        if ier:
            msg = "An error occurred after calling MINTLIB.mnt_extensivefieldconverter_getFaceData."
            warning_handler(FILE, 'getFaceData', ier,
                            detailedmsg=msg)
        return res

    def getCellByCellDataFromUniqueEdgeData(self, unique_edge_data):
        """
        Get the cell by cell data from the unique edge data
        :param unique_edge_data: data defined on unique dege Ids, size of array should be numEdges
        :returns cell by cell data
        """
        n = self.numGridCells * NUM_EDGES_PER_QUAD
        res = numpy.empty(n, numpy.float64)
        MINTLIB.mnt_extensivefieldconverter_getCellByCellDataFromUniqueEdgeData.argtypes = [
                                                                                POINTER(c_void_p),
                                                                                DOUBLE_ARRAY_PTR,
                                                                                DOUBLE_ARRAY_PTR]
        ier = MINTLIB.mnt_extensivefieldconverter_getCellByCellDataFromUniqueEdgeData(self.obj, unique_edge_data, res)
        if ier:
            msg = "An error occurred after calling MINTLIB.mnt_extensivefieldconverter_getCellByCellDataFromUniqueEdgeData."
            warning_handler(FILE, 'getCellByCellDataFromUniqueEdgeData', ier,
                            detailedmsg=msg)
        return res

