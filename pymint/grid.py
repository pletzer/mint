from ctypes import c_void_p, c_double, c_int, byref, POINTER, c_char_p, c_size_t, c_longlong
from . import LIB
import numpy

def error_handler(filename, methodname, ier):
    raise RuntimeError(f'ERROR ier={ier} after calling {methodname} in {filename}!')

FILE = 'grid.py'

class Grid(object):


    def __init__(self):
        """
        Regrid edge field constructor.
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        LIB.mnt_grid_new.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_grid_new(self.obj)
        if ier: 
            error_handler(FILE, '__init__', ier)


    def __del__(self):
        """
        Regrid edge field destructor.
        """
        LIB.mnt_grid_del.argtypes = [POINTER(c_void_p)]
        ier = LIB.mnt_grid_del(self.obj)
        if ier: 
            error_handler(FILE, '__del__', ier)


    def setFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length should be added/subtracted
                                     in order to make eeach cell as compact as poassible
        :param averageLonAtPole: set to 1 if the longitudes at the poles should the average of 
                                 the cell's longitudes

        note:: a lon-lat grid requires 0, 0 and a cibed sphere grid requires 1, 1
        """
        LIB.mnt_grid_setFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = LIB.mnt_grid_setFlags(self.obj, fixLonAcrossDateline, averageLonAtPole)
        if ier:
            error_handler(FILE, 'setFlags', ier)


    def loadFrom2DUgrid(self, fileAndMeshName):
        """
        Load a grid from a 2D UGRID file.

        :param fileAndMeshName: string in the format filename:meshname
        """
        LIB.mnt_grid_loadFrom2DUgrid.argtypes = [POINTER(c_void_p), c_char_p]
        fm = fileAndMeshName.encode('utf-8')
        ier = LIB.mnt_grid_loadFrom2DUgrid(self.obj, fm)
        if ier:
            error_handler(FILE, 'loadFrom2DUgrid', ier)


    def load(self, filename):
        """
        Load the grid from a VTK file.

        :param filename: file name
        """
        LIB.mnt_grid_load.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = LIB.mnt_grid_load(self.obj, fm)
        if ier:
            error_handler(FILE, 'load', ier)


    def dump(self, filename):
        """
        Dump the grid to a VTK file.

        :param filename: file name
        """
        LIB.mnt_grid_dump.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = LIB.mnt_grid_dump(self.obj, fm)
        if ier:
            error_handler(FILE, 'dump', ier)



    def setPoints(self, points):
        """
        Set the points coordinates

        :param points: numpy contiguous array of shape (ncells, num_verts_per_cell, 3)
        """
        ncells, num_verts_per_cell, ndim = points.shape
        if ndim != 3:
            raise RuntimeError(f'ERROR: points.shape[2] should be 3, got {ndim}!')
        LIB.mnt_grid_setPointsPtr.argtypes = [POINTER(c_void_p), c_int, c_longlong, 
                                              numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        ier = LIB.mnt_grid_setPointsPtr(self.obj, num_verts_per_cell, ncells, points)
        if ier:
            error_handler(FILE, 'setPointsPtr', ier)

    def getEdgeId(self, cellId, edgeIndex):
        """
        Get the edge Id and direction of a cellId, edgeIndex pair

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns an edge index, sign pair
        """
        LIB.mnt_grid_getEdgeId.argtypes = [POINTER(c_void_p), c_longlong, c_int, POINTER(c_size_t), POINTER(c_int)]
        edgeId = c_size_t()
        edgeSign = c_int()
        ier = LIB.mnt_grid_getEdgeId(self.obj, cellId, edgeIndex, byref(edgeId), byref(edgeSign))
        if ier:
            error_handler(FILE, 'getEdgeId', ier)
        return (edgeId.value, edgeSign.value)


    def getNodeIds(self, cellId, edgeIndex):
        """
        Get the node Ids of a cellId, edgeIndex pair

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns two node indices
        """
        LIB.mnt_grid_getNodeIds.argtypes = [POINTER(c_void_p), c_longlong, c_int, POINTER(c_size_t)]
        nodeIds = (c_size_t*2)()
        ier = LIB.mnt_grid_getNodeIds(self.obj, cellId, edgeIndex, nodeIds)
        if ier:
            error_handler(FILE, 'getNodeIds', ier)
        return (nodeIds[0], nodeIds[1])


    def attach(self, varname, data):
        """
        Attach data to the grid.

        :param varname: field name
        :param data: numpy array of size (ncells, num_edges_per_cell, ndims)
        """
        ndims = data.shape[-1]
        ncells = self.getNumberOfCells()
        nEdgesPerCell = 4 # self.getNumEdgesPerCell()
        nDataPerCell = nEdgesPerCell * ndims
        LIB.mnt_grid_attach.argtypes = [POINTER(c_void_p), c_char_p, c_int, 
                                        numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        ier = LIB.mnt_grid_attach(self.obj, varname.encode('utf-8'), nDataPerCell, data)
        if ier:
            error_handler(FILE, 'attach', ier)
            

    def getNumberOfCells(self):
        """
        Get the number of cells

        :returns number
        """
        LIB.mnt_grid_getNumberOfCells.argtypes = [POINTER(c_void_p), POINTER(c_size_t)]
        n = c_size_t()
        ier = LIB.mnt_grid_getNumberOfCells(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumberOfCells', ier)
        return n.value


    def getNumberOfEdges(self):
        """
        Get the number of unique edges of the grid.

        :returns number
        """
        LIB.mnt_grid_getNumberOfEdges.argtypes = [POINTER(c_void_p)]
        n = c_size_t()
        ier = LIB.mnt_grid_getNumberOfEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumberOfEdges', ier)
        return n.value

