from ctypes import (c_void_p, c_int, byref, POINTER, c_char_p,
                    c_size_t, c_longlong, c_double)
from . import MINTLIB, NUM_VERTS_PER_QUAD, NUM_VERTS_PER_EDGE
from . import error_handler
import numpy


FILE = 'grid.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)
SIZE_T_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.uintp)


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

    def setFlags(self, fixLonAcrossDateline, averageLonAtPole, degrees=True) :
        """
        Set the grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length should be
                                     added/subtracted to make each cell as
                                     compact as possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles should
                                 be the average of the cell's longitudes
        :param degrees: set to False if coordinates are note in degrees

        note:: a lon-lat grid requires fixLonAcrossDateline=0, averageLonAtPole=0 and 
               cubed sphere grid requires fixLonAcrossDateline=1, averageLonAtPole=1
        note:: call this before reading the grid from file
        """
        MINTLIB.mnt_grid_setFlags.argtypes = [POINTER(c_void_p), c_int, c_int, c_int]
        degrees_int = int(degrees)
        ier = MINTLIB.mnt_grid_setFlags(self.obj,
                                    fixLonAcrossDateline,
                                    averageLonAtPole,
                                    degrees_int)
        if ier:
            error_handler(FILE, 'setFlags', ier)

    def loadFromUgrid2DData(self, xyz: numpy.ndarray, face2nodes: numpy.ndarray, edge2nodes: numpy.ndarray) -> None:
        """
        Load a grid from 2D UGRID data structures.

        :param xyz: array of vertex coordinates, size npoints * 3
        :param face2nodes: array of face/cell connectivity to vertex Ids
        :param edge2nodes: array of edge connectivity to vertex Ids
        """
        if not xyz.flags.c_contiguous or not face2nodes.flags.c_contiguous or not edge2nodes.flags.c_contiguous:
            error_handler(FILE, 'loadFromUgrid2DData', 'input arrays must be C contiguous!')

        MINTLIB.mnt_grid_loadFromUgrid2DData.argtypes = [POINTER(c_void_p),
                                                         c_size_t, c_size_t, c_size_t,
                                                         DOUBLE_ARRAY_PTR,
                                                         SIZE_T_ARRAY_PTR, SIZE_T_ARRAY_PTR]
        ncells = c_size_t(face2nodes.shape[0])
        nedges = c_size_t(edge2nodes.shape[0])
        npoints = c_size_t(xyz.shape[0])
        ier = MINTLIB.mnt_grid_loadFromUgrid2DData(self.obj, ncells, nedges, npoints,
                                                   xyz, face2nodes, edge2nodes)
        if ier:
            error_handler(FILE, 'loadFromUgrid2DData', ier)

    def loadFromUgrid2DFile(self, fileAndMeshName: str) -> None:
        """
        Load a grid from a 2D UGRID file.

        :param fileAndMeshName: string in the format filename$meshname
        """
        MINTLIB.mnt_grid_loadFromUgrid2DFile.argtypes = [POINTER(c_void_p), c_char_p]
        fm = fileAndMeshName.encode('utf-8')
        ier = MINTLIB.mnt_grid_loadFromUgrid2DFile(self.obj, fm)
        if ier:
            error_handler(FILE, 'loadFromUgrid2DFile', ier)

    def loadFromFV32DData(self, lon_corners: numpy.ndarray, lat_corners: numpy.ndarray) -> None:
        """
        Load a grid from FV3/SCRIP-style cell-corner arrays.

        :param lon_corners: longitudes of cell corners, shape (ncells, 4), degrees, C contiguous
        :param lat_corners: latitudes  of cell corners, shape (ncells, 4), degrees, C contiguous

        .. note:: corner ordering is SW, SE, NE, NW (CCW)
        .. note:: no edge or face-node connectivity is built; getNodeIds/getEdgeId will fail
        .. note:: call setFlags before this method to control dateline and pole handling
        """
        if not lon_corners.flags.c_contiguous or not lat_corners.flags.c_contiguous:
            error_handler(FILE, 'loadFromFV32DData', 'input arrays must be C contiguous!')

        ncells = lon_corners.shape[0]
        if lon_corners.shape != (ncells, 4) or lat_corners.shape != (ncells, 4):
            error_handler(FILE, 'loadFromFV32DData',
                          f'expected shape (ncells, 4), got {lon_corners.shape} / {lat_corners.shape}')

        MINTLIB.mnt_grid_loadFromFV32DData.argtypes = [POINTER(c_void_p),
                                                       c_size_t,
                                                       DOUBLE_ARRAY_PTR,
                                                       DOUBLE_ARRAY_PTR]
        ier = MINTLIB.mnt_grid_loadFromFV32DData(self.obj, c_size_t(ncells),
                                                  lon_corners, lat_corners)
        if ier:
            error_handler(FILE, 'loadFromFV32DData', ier)

    def loadFromFV32DFile(self, filename: str) -> None:
        """
        Load a grid from a FV3/SCRIP NetCDF file.

        :param filename: path to a SCRIP-format NetCDF file containing
                         ``grid_corner_lon`` and ``grid_corner_lat`` variables
                         of dimensions ``(grid_size, grid_corners=4)``

        .. note:: call setFlags before this method to control dateline and pole handling
        """
        MINTLIB.mnt_grid_loadFromFV32DFile.argtypes = [POINTER(c_void_p), c_char_p]
        ier = MINTLIB.mnt_grid_loadFromFV32DFile(self.obj, filename.encode('utf-8'))
        if ier:
            error_handler(FILE, 'loadFromFV32DFile', ier)

    def load(self, filename: str):
        """
        Load the grid from a VTK file.

        :param filename: file name
        """
        MINTLIB.mnt_grid_load.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = MINTLIB.mnt_grid_load(self.obj, fm)
        if ier:
            error_handler(FILE, 'load', ier)

    def dump(self, filename: str) -> None:
        """
        Dump the grid to a VTK file.

        :param filename: file name
        """
        MINTLIB.mnt_grid_dump.argtypes = [POINTER(c_void_p), c_char_p]
        fm = filename.encode('utf-8')
        ier = MINTLIB.mnt_grid_dump(self.obj, fm)
        if ier:
            error_handler(FILE, 'dump', ier)

    def setPoints(self, points: numpy.ndarray) -> None:
        """
        Set the points (vertices), cell by cell. The points should ordered in counterclockwise way.

        :param points: numpy contiguous array of shape (ncells,
                       num_verts_per_cell, 3) using C ordering

        note:: mnt_grid_setPointsPtr hands the C++ side a raw pointer into
               `points`, it does not copy the data (same as attach(...,
               copy=False)). We keep our own reference to `points` here so
               it stays alive for as long as this Grid does -- otherwise,
               if the caller lets their own reference go out of scope (e.g.
               a helper function that builds the array and returns just the
               Grid), Python is free to garbage-collect it, leaving the C++
               side holding a dangling pointer into freed/reused memory.
        """
        if not points.flags.c_contiguous:
            error_handler(FILE, 'setPoints', 'points array must be C contiguous!')

        ncells, num_verts_per_cell, ndim = points.shape
        if ndim != 3:
            error_handler(FILE, f'setPoints: points.shape[2] != 3, got {ndim}!')
        MINTLIB.mnt_grid_setPointsPtr.argtypes = [POINTER(c_void_p),
                                              DOUBLE_ARRAY_PTR]
        MINTLIB.mnt_grid_build.argtypes = [POINTER(c_void_p), c_int, c_longlong]
        ier = MINTLIB.mnt_grid_setPointsPtr(self.obj, points)
        if ier:
            error_handler(FILE, 'setPointsPtr', ier)
        ier = MINTLIB.mnt_grid_build(self.obj, num_verts_per_cell, ncells)
        if ier:
            error_handler(FILE, 'build', ier)
        # keep a reference to the points array so it is not garbage collected
        self._points = points

    def getPoints(self) -> numpy.ndarray:
        """
        Get a view of the point (vertex) array of the cell-by-cell mesh.

        :return numpy contiguous array of shape (ncells,
                      num_verts_per_cell, 3)

        note: the returned array is a view into the C++ side's memory, so it will be 
        invalidated if the C++ side is deleted or if setPoints is called again. If
        """
        MINTLIB.mnt_grid_getPointsPtr.argtypes = [POINTER(c_void_p), POINTER(POINTER(c_double))]

        pointsPtr = POINTER(c_double)()
        ier = MINTLIB.mnt_grid_getPointsPtr(self.obj, byref(pointsPtr))
        if ier:
            error_handler(FILE, 'getPointsPtr', ier)
            
        ncells = self.getNumberOfCells()
        # create a numpy array from a C pointer
        return numpy.ctypeslib.as_array(pointsPtr, shape=(ncells, NUM_VERTS_PER_QUAD, 3))

    def getEdgeId(self, cellId: int, edgeIndex: int) -> tuple[int, int]:
        """
        Get the edge Id and direction of a cellId, edgeIndex pair.

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns a unique edge Id, sign pair
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

    def getNodeIds(self, cellId: int, edgeIndex: int) -> tuple[int, int]:
        """
        Get the node Ids of a cellId, edgeIndex pair.

        :param cellId: Id of the cell
        :param edgeIndex: edge index of the cell (0...3)
        :returns node indices
        """
        MINTLIB.mnt_grid_getNodeIds.argtypes = [POINTER(c_void_p),
                                                c_size_t, c_int,
                                                POINTER(c_size_t)]
        nodeIds = (c_size_t*NUM_VERTS_PER_EDGE)()
        ier = MINTLIB.mnt_grid_getNodeIds(self.obj, cellId, edgeIndex, nodeIds)
        if ier:
            error_handler(FILE, 'getNodeIds', ier)
        return nodeIds[0], nodeIds[1]

    def computeEdgeArcLengths(self) -> None:
        """
        Compute and store edge arc lengths.
        :note assumes the sphere radius to be one
        """
        MINTLIB.mnt_grid_computeEdgeArcLengths.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_grid_computeEdgeArcLengths(self.obj)
        if ier:
            error_handler(FILE, 'computeEdgeArcLengths', ier)

    def getEdgeArcLength(self, cellId: int, edgeIndex: int) -> float:
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

    def attach(self, varname: str, data: numpy.ndarray, copy: bool=True):
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

    def getNumberOfCells(self) -> int:
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

    def getNumberOfEdges(self) -> int:
        """
        Get the number of unique edges of the grid.

        :returns number
        """
        MINTLIB.mnt_grid_getNumberOfEdges.argtypes = [POINTER(c_void_p),
                                                      POINTER(c_size_t)]
        n = c_size_t(0)
        ier = MINTLIB.mnt_grid_getNumberOfEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumberOfEdges', ier)
        return n.value

    def check(self) -> int:
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


