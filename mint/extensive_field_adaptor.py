from ctypes import (c_void_p, c_int, byref, POINTER)
from . import MINTLIB, CELL_BY_CELL_DATA, NUM_EDGES_PER_QUAD, FUNC_SPACE_W1, FUNC_SPACE_W2
from . import error_handler, warning_handler
import numpy


FILE = 'extensive_field_adaptor.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


class ExtensiveFieldAdaptor(object):
    """
    A class to convert extensive fields to vector fields and vice versa
    """

    def __init__(self):
        """
        Constructor.
        """

        self.numGridEdges = 0
        self.numGridCells = 0
        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        MINTLIB.mnt_extensivefieldadaptor_new.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_extensivefieldadaptor_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_extensivefieldadaptor_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_extensivefieldadaptor_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setGrid(self, grid):
        """
        Set the grid.

        :param grid: instance of Grid
        """
        self.numGridEdges = grid.getNumberOfEdges()
        self.numGridCells = grid.getNumberOfCells()
        MINTLIB.mnt_extensivefieldadaptor_setGrid.argtypes = [POINTER(c_void_p), c_void_p]
        ier = MINTLIB.mnt_extensivefieldadaptor_setGrid(self.obj, grid.ptr)
        if ier:
            error_handler(FILE, 'setGrid', ier)

    def fromVectorField(self, u, v, data, placement, fs):
        """
        Get the edge integrated data from a vector field
        :param u: zonal (x) component of vectors on edges, array (see placement argument below)
        :param v: meridonal (y) component of vectors on edges, array (see placement argument below)
        :param data: extensive field (output), array (see placement argument below)
        :param placement: mint.CELL_BY_CELL_DATA if vx and vy are cell by cell (size num cells * mint.NUM_EDGES_PER_QUAD),
                          vx and vy are assumed to be on unique edges otherwise (size num edges)
        :param fs: function space, either FUNC_SPACE_W1 or FUNC_SPACE_W2
        :note: A radius of 1 is assumed for the planet.
        """

        # check the function space
        if fs not in (FUNC_SPACE_W1, FUNC_SPACE_W2):
             msg = f"function space must be either FUNC_SPACE_W1 or FUNC_SPACE_W2 (!= {fs})"
             ier = 11
             error_handler(FILE, 'fromVectorField', ier, detailedmsg=msg)

        # check the data size
        n0 = self.numGridCells * NUM_EDGES_PER_QUAD
        if placement != CELL_BY_CELL_DATA:
            n0 = self.numGridEdges

        if u.shape[-1] != n0 or v.shape[-1] != n0 or data.shape[-1] != n0:
            msg = f"u, v have wrong size (= {u.shape[-1]}, {v.shape[-1], {data.shape[-1]}} != {n0})"
            ier = 10
            error_handler(FILE, 'fromVectorField', ier, detailedmsg=msg)

        res = numpy.empty(n0, numpy.float64)

        MINTLIB.mnt_extensivefieldadaptor_fromVectorField.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int,
                                                        c_int]
        ier = MINTLIB.mnt_extensivefieldadaptor_fromVectorField(self.obj, u, v, data, placement, fs)
        if ier:
            msg = "An error occurred after calling MINTLIB.mnt_extensivefieldadaptor_fromVectorField."
            warning_handler(FILE, 'fromVectorField', ier,
                            detailedmsg=msg)


    def toVectorField(self, edge_data, face_data, u, v, placement):
        """
        Compute the vector components from the edge and face extensive fields
        :param edge_data: extensive edge data
        :param face_data: extensive face data
        :param u: zonal (x) component (output)
        :param v: meridional (y) component (output)
        :param placement: mint.CELL_BY_CELL_DATA if the data are cell by cell (size num cells * mint.NUM_EDGES_PER_QUAD),
                          otherwise assumed to be on unique edges (size num edges)
        :note: a radius of 1 is assumed for the planet.
        """

        # check the data size
        n0 = self.numGridCells * NUM_EDGES_PER_QUAD
        if placement != CELL_BY_CELL_DATA:
            n0 = self.numGridEdges

        if edge_data.shape[-1] != n0 or face_data.shape[-1] != n0 or u.shape[-1] != n0 or v.shape[-1] != n0:
            msg = f"""supplied edge or face data have wrong size
(= {edge_data.shape[-1]}, {face_data.shape[-1]}, {u.shape[-1]}, {v.shape[-1]} != {n0})"""
            ier = 10
            error_handler(FILE, 'toVectorField', ier, detailedmsg=msg)

        MINTLIB.mnt_extensivefieldadaptor_toVectorField.argtypes = [
                                                        POINTER(c_void_p),
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        DOUBLE_ARRAY_PTR,
                                                        c_int]
        ier = MINTLIB.mnt_extensivefieldadaptor_toVectorField(self.obj, edge_data, face_data, u, v, placement)
        if ier:
            msg = "An error occurred after calling MINTLIB.mnt_extensivefieldadaptor_toVectorField."
            warning_handler(FILE, 'toVectorField', ier,
                            detailedmsg=msg)

