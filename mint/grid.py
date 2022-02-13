from ctypes import (c_void_p, c_int, byref, POINTER, c_char_p,
                    c_size_t, c_longlong, c_double)
from . import MINTLIB
from . import error_handler, warning_handler
import numpy


FILE = 'grid.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class Grid(object):
    """
    A class to represent a collection of quad cells.
    """

    def __init__(self):
        """
        Constructor.
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)
        # container holding attached data
        self.data = {}

        MINTLIB.mnt_grid_new.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_grid_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_grid_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_grid_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length should be
                                     added/subtracted to make each cell as
                                     compact as possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles should
                                 be the average of the cell's longitudes

        note:: a lon-lat grid requires 0, 0 and cubed sphere grid requires 1, 1
        note:: call this before reading the grid from file
        """
        MINTLIB.mnt_grid_setFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = MINTLIB.mnt_grid_setFlags(self.obj,
                                    fixLonAcrossDateline,
                                    averageLonAtPole)
        if ier:
            error_handler(FILE, 'setFlags', ier)

    def loadFromUgrid2D(self, fileAndMeshName):
        """
        Load a grid from a 2D UGRID file.

        :param fileAndMeshName: string in the format filename$meshname
        """
        MINTLIB.mnt_grid_loadFromUgrid2D.argtypes = [POINTER(c_void_p), c_char_p]
        fm = fileAndMeshName.encode('utf-8')
        ier = MINTLIB.mnt_grid_loadFromUgrid2D(self.obj, fm)
        if ier:
            error_handler(FILE, 'loadFromUgrid2D', ier)

    def load(self, filename):
        """
        Load the grid from a VTK file.

        :param filename: file name
        """
        MINTLIB.mnt_grid_load.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = MINTLIB.mnt_grid_load(self.obj, fm)
        if ier:
            error_handler(FILE, 'load', ier)

    def dump(self, filename):
        """
        Dump the grid to a VTK file.

        :param filename: file name
        """
        MINTLIB.mnt_grid_dump.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = MINTLIB.mnt_grid_dump(self.obj, fm)
        if ier:
            error_handler(FILE, 'dump', ier)

    def setPoints(self, points):
        """
        Set the point (vertices) and build the connectivity.

        :param points: numpy contiguous array of shape (ncells,
                       num_verts_per_cell, 3) using C ordering
        """
        ncells, num_verts_per_cell, ndim = points.shape
        if ndim != 3:
            raise RuntimeError(f'ERROR: points.shape[2] != 3, got {ndim}!')
        MINTLIB.mnt_grid_setPointsPtr.argtypes = [POINTER(c_void_p),
                                              DOUBLE_ARRAY_PTR]
        MINTLIB.mnt_grid_build.argtypes = [POINTER(c_void_p), c_int, c_longlong]
        ier = MINTLIB.mnt_grid_setPointsPtr(self.obj, points)
        if ier:
            error_handler(FILE, 'setPointsPtr', ier)
        ier = MINTLIB.mnt_grid_build(self.obj, num_verts_per_cell, ncells)
        if ier:
            error_handler(FILE, 'setPointsPtr', ier)

    def getPoints(self):
        """
        Get a view of the points (vertices) array of the cell-by-cell mesh.

        :return numpy contiguous array of shape (ncells,
                      num_verts_per_cell, 3)
        """
        MINTLIB.mnt_grid_getPointsPtr.argtypes = [POINTER(c_void_p), POINTER(POINTER(c_double))]

        pointsPtr = POINTER(c_double)()
        ier = MINTLIB.mnt_grid_getPointsPtr(self.obj, byref(pointsPtr))
        if ier:
            error_handler(FILE, 'getPointsPtr', ier)
        # create a numpy array from a C pointer\
        ncells = self.getNumberOfCells()
        return numpy.ctypeslib.as_array(pointsPtr, shape=(ncells, 4, 3))

    def getEdgeId(self, cellId, edgeIndex):
        """
        Get the edge Id and direction of a cellId, edgeIndex pair.

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns an edge index, sign pair
        """
        MINTLIB.mnt_grid_getEdgeId.argtypes = [POINTER(c_void_p),
                                           c_longlong, c_int,
                                           POINTER(c_size_t),
                                           POINTER(c_int)]
        edgeId = c_size_t()
        edgeSign = c_int()
        ier = MINTLIB.mnt_grid_getEdgeId(self.obj, cellId, edgeIndex,
                                     byref(edgeId), byref(edgeSign))
        if ier:
            error_handler(FILE, 'getEdgeId', ier)
        return edgeId.value, edgeSign.value

    def getNodeIds(self, cellId, edgeIndex):
        """
        Get the node Ids of a cellId, edgeIndex pair.

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns two node indices
        """
        MINTLIB.mnt_grid_getNodeIds.argtypes = [POINTER(c_void_p),
                                            c_longlong, c_int,
                                            POINTER(c_size_t)]
        nodeIds = (c_size_t*2)()
        ier = MINTLIB.mnt_grid_getNodeIds(self.obj, cellId, edgeIndex, nodeIds)
        if ier:
            error_handler(FILE, 'getNodeIds', ier)
        return nodeIds[0], nodeIds[1]

    def computeEdgeArcLengths(self):
        """
        Compute and store edge arc lengths.
        :note assumes the sphere radius to be one
        """
        MINTLIB.mnt_grid_computeEdgeArcLengths.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_grid_computeEdgeArcLengths(self.obj)
        if ier:
            error_handler(FILE, 'computeEdgeArcLengths', ier)

    def getEdgeArcLength(self, cellId, edgeIndex):
        """
        Get the arc length for given cell and edge.
        :param cellId: cell Id
        :param edgeIndex: edge index (0...3)
        :returns length assuming radius of one
        """
        res = c_double()
        MINTLIB.mnt_grid_getEdgeArcLength.argtypes = [POINTER(c_void_p),
                                                  c_longlong, c_int,
                                                  POINTER(c_double)]
        ier = MINTLIB.mnt_grid_getEdgeArcLength(self.obj, cellId, edgeIndex,
                                            byref(res))
        if ier:
            error_handler(FILE, 'getEdgeArcLength', ier)
        return res.value

    def attach(self, varname, data, copy=True):
        """
        Attach data to the grid.

        :param varname: field name
        :param data: numpy array of size (ncells, :)
        :param copy: set to False if you do not want the data to be copied
        """
        nDataPerCell = 1
        if len(data.shape) > 1:
            nDataPerCell = data.shape[-1]
        MINTLIB.mnt_grid_attach.argtypes = [POINTER(c_void_p), c_char_p, c_int,
                                            DOUBLE_ARRAY_PTR]
        if copy:
            # make a copy to ensure that the data exist during the life of
            # this instance
            self.data[varname] = data.copy()
            ier = MINTLIB.mnt_grid_attach(self.obj, varname.encode('utf-8'),
                                          nDataPerCell, self.data[varname])
        else:
            ier = MINTLIB.mnt_grid_attach(self.obj, varname.encode('utf-8'),
                                          nDataPerCell, data)
        if ier:
            error_handler(FILE, 'attach', ier)

    def getNumberOfCells(self):
        """
        Get the number of cells.

        :returns number
        """
        MINTLIB.mnt_grid_getNumberOfCells.argtypes = [POINTER(c_void_p),
                                                      POINTER(c_size_t)]
        n = c_size_t()
        ier = MINTLIB.mnt_grid_getNumberOfCells(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumberOfCells', ier)
        return n.value

    def getNumberOfEdges(self):
        """
        Get the number of unique edges of the grid.

        :returns number
        """
        MINTLIB.mnt_grid_getNumberOfEdges.argtypes = [POINTER(c_void_p),
                                                      POINTER(c_size_t)]
        n = c_size_t()
        ier = MINTLIB.mnt_grid_getNumberOfEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumberOfEdges', ier)
        elif n.value <= 0:
            msg = f'''A number of unique edges ({n}) is expected if the grid was not loaded from UGRID,
for instance if you set the grid points cell by cell using the setPoints method'''
            warning_handler(FILE, 'getNumberOfEdges', ier, detailedmsg=msg)
        return n.value

    def check(self):
        """
        Check that the cells have positive area.

        :returns number of bad cells
        """
        num_bad_cells = c_size_t()
        MINTLIB.mnt_grid_check.argtypes = [POINTER(c_void_p), POINTER(c_size_t)]
        ier = MINTLIB.mnt_grid_check(self.obj, byref(num_bad_cells))
        if ier:
            error_handler(FILE, 'check', ier)
        return num_bad_cells.value


