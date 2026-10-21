from ctypes import (c_void_p, c_int, c_size_t, byref, POINTER)
from . import MINTLIB, CELL_BY_CELL_DATA
from . import error_handler
import numpy


FILE = 'whitney_element.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64, flags="C_CONTIGUOUS")
LONGLONG_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.int64, flags="C_CONTIGUOUS")


class WhitneyElement(object):
    """
    An edge (1-form) element evaluated directly, outside of VectorInterp,
    PolylineIntegral and RegridEdges. A field represented by its edge
    integrals is returned as covariant components (a, b), i.e. the 1-form
    a dxi + b deta in the parametric coordinates (xi, eta) of a grid cell.
    """

    def __init__(self, order, grid, biquadratic=True):
        """
        Constructor.

        :param order: 1 for the lowest order Whitney element, 2 for the second order
                      element built on 2x2 blocks of grid cells (same edge data)
        :param grid: instance of Grid, the grid must have UGRID connectivity for order 2
        :param biquadratic: order 2 only. True to map each 2x2 block with its
                            biquadratic 9-node map, False to use the piecewise bilinear
                            cell coordinates X = (p + xi)/2, Y = (q + eta)/2, which are
                            exact when the field is reconstructed in coordinates in which
                            the grid is uniform (e.g. the panel chart of an equiangular
                            cubed sphere)
        """
        self.ptr = c_void_p()
        self.obj = byref(self.ptr)
        self.numGridCells = grid.getNumberOfCells()
        MINTLIB.mnt_whitneyelement_new.argtypes = [POINTER(c_void_p), c_int, c_void_p, c_int]
        ier = MINTLIB.mnt_whitneyelement_new(self.obj, order, grid.ptr, int(biquadratic))
        if ier:
            error_handler(FILE, '__init__', ier)
        # the C side borrows the grid
        self._grid = grid

    def __del__(self):
        """
        Destructor.
        """
        if self.ptr:
            MINTLIB.mnt_whitneyelement_del.argtypes = [POINTER(c_void_p)]
            MINTLIB.mnt_whitneyelement_del(self.obj)

    def getOrder(self):
        """
        Get the order of the element.

        :returns: 1 or 2
        """
        MINTLIB.mnt_whitneyelement_getOrder.argtypes = [POINTER(c_void_p), POINTER(c_int)]
        order = c_int()
        ier = MINTLIB.mnt_whitneyelement_getOrder(self.obj, byref(order))
        if ier:
            error_handler(FILE, 'getOrder', ier)
        return order.value

    def getCovariantComponents(self, cellIds, xi, eta, data, placement=CELL_BY_CELL_DATA):
        """
        Get the covariant components of a field at points given by their cell and
        parametric coordinates.

        :param cellIds: cell Id of each point
        :param xi: first parametric coordinate of each point in its cell
        :param eta: second parametric coordinate of each point in its cell
        :param data: edge integrals, cell by cell (numCells*4, each cell's edges in the
                     positive xi/eta direction) or on unique edges
        :param placement: CELL_BY_CELL_DATA or UNIQUE_EDGE_DATA
        :returns: arrays a, b of the 1-form a dxi + b deta at the points
        :note: an error is raised if any cell is invalid
        """
        cellIds = numpy.ascontiguousarray(cellIds, numpy.int64)
        xi = numpy.ascontiguousarray(xi, numpy.float64)
        eta = numpy.ascontiguousarray(eta, numpy.float64)
        data = numpy.ascontiguousarray(data, numpy.float64)
        n = cellIds.size
        if xi.size != n or eta.size != n:
            error_handler(FILE, 'getCovariantComponents', -1,
                          'cellIds, xi and eta must have the same size')
        a = numpy.empty(n, numpy.float64)
        b = numpy.empty(n, numpy.float64)
        MINTLIB.mnt_whitneyelement_getCovariantComponents.argtypes = [
            POINTER(c_void_p), c_void_p, c_size_t, LONGLONG_ARRAY_PTR, DOUBLE_ARRAY_PTR,
            DOUBLE_ARRAY_PTR, DOUBLE_ARRAY_PTR, c_int, DOUBLE_ARRAY_PTR, DOUBLE_ARRAY_PTR]
        ier = MINTLIB.mnt_whitneyelement_getCovariantComponents(self.obj, self._grid.ptr, n,
                                                                cellIds, xi, eta, data,
                                                                placement, a, b)
        if ier:
            error_handler(FILE, 'getCovariantComponents', ier,
                          f'{ier} points have an invalid cell')
        return a.reshape(numpy.shape(cellIds)), b.reshape(numpy.shape(cellIds))
