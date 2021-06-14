from ctypes import c_void_p, c_double, c_int, byref, POINTER,\
                   c_char_p, c_size_t
from . import LIB
import numpy
import logging


FILE = 'regrid_edges.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)


def error_handler(filename, methodname, ier):
    msg = f'ier={ier} after calling {methodname} in {filename}!'
    logging.error(msg)
    raise RuntimeError(msg)


def warning_handler(filename, methodname, ier, detailedMsg=''):
    msg = f'ier={ier} after calling {methodname} in {filename}!\n'
    msg += detailedMsg
    logging.warning(msg)


class RegridEdges(object):
    """
    A class to interpolate a 1 form from a grid to another grid
    """

    def __init__(self):
        """
        Regrid edge field constructor.
        """

        self.obj = byref(c_void_p())

        ier = LIB.mnt_regridedges_new(self.obj)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Regrid edge field destructor.
        """
        ier = LIB.mnt_regridedges_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def setSrcGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the source grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length
                                     should be added/subtracted
                                     to make each cell as compact as
                                     possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles
                                 should be the average of the cells'
                                 longitudes

        note:: a lon-lat grid requires 0, 0 and a cubed sphere grid
               requires 1, 1
        """
        LIB.mnt_regridedges_setSrcGridFlags.argtypes = [POINTER(c_void_p),
                                                        c_int, c_int]
        ier = LIB.mnt_regridedges_setSrcGridFlags(self.obj,
                                                  fixLonAcrossDateline,
                                                  averageLonAtPole)
        if ier:
            error_handler(FILE, 'setSrcGridFlags', ier)

    def setDstGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the destrination grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length
                                     should be added/subtracted
                                     to make each cell as compact as
                                     possible
        :param averageLonAtPole: set to 1 if the longitudes at the poles
                                 should be the average of the cells'
                                 longitudes

        note:: a lon-lat grid requires 0, 0 and a cubed sphere grid
               requires 1, 1
        """
        LIB.mnt_regridedges_setDstGridFlags.argtypes = [POINTER(c_void_p),
                                                        c_int, c_int]
        ier = LIB.mnt_regridedges_setDstGridFlags(self.obj,
                                                  fixLonAcrossDateline,
                                                  averageLonAtPole)
        if ier:
            error_handler(FILE, 'setDstGridFlags', ier)

    def loadSrcGrid(self, filename):
        """
        Load a source grid from a 2D UGRID file.

        :param filename: string in the format filename:meshname
        """
        LIB.mnt_regridedges_loadSrcGrid.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = LIB.mnt_regridedges_loadSrcGrid(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'loadSrcGrid', ier)

    def loadDstGrid(self, filename):
        """
        Load a destination grid from a 2D UGRID file.

        :param filename: string in the format filename:meshname
        """
        LIB.mnt_regridedges_loadDstGrid.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = LIB.mnt_regridedges_loadDstGrid(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'loadDstGrid', ier)

    def getNumSrcEdges(self):
        """
        Get the number of unique edges of the source grid.

        :returns number
        """
        LIB.mnt_regridedges_getNumSrcEdges.argtypes = [POINTER(c_void_p)]
        LIB.mnt_regridedges_getNumSrcEdges.restype = c_size_t
        n = c_size_t()
        ier = LIB.mnt_regridedges_getNumSrcEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumSrcEdges', ier)
        return n.value

    def getNumDstEdges(self):
        """
        Get the number of unique edges of the destination grid.

        :returns number
        """
        LIB.mnt_regridedges_getNumDstEdges.argtypes = [POINTER(c_void_p)]
        LIB.mnt_regridedges_getNumDstEdges.restype = c_size_t
        n = c_size_t()
        ier = LIB.mnt_regridedges_getNumDstEdges(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumDstEdges', ier)
        return n.value

    def build(self, numCellsPerBucket, periodX, debug):
        """
        Build the regridder and compute the regridding weights

        :param numCellsPerBucket: average number of cells per bucket
                                  (performance only)
        :param periodX: periodicity length (set to 0 if non-periodic)
        :param debug: 0=no debug info, 1=print debug info, 2=save bad
                      edges in VTK file
        """
        LIB.mnt_regridedges_build.argtypes = [POINTER(c_void_p), c_int,
                                              c_double, c_int]
        ier = LIB.mnt_regridedges_build(self.obj, numCellsPerBucket,
                                        periodX, debug)
        if ier:
            msg = "Some target lines fall outside the grid."
            warning_handler(FILE, 'build', ier, detailedMsg=msg)

    def dumpWeights(self, filename):
        """
        Dump the weights to a file
        :param filename: file name
        """
        LIB.mnt_regridedges_dumpWeights.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = LIB.mnt_regridedges_dumpWeights(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'dumpWeights', ier)

    def loadWeights(self, filename):
        """
        Load the weights from a file
        :param filename: file name
        """
        LIB.mnt_regridedges_loadWeights.argtypes = [POINTER(c_void_p),
                                                    c_char_p, c_int]
        fn = filename.encode('utf-8')
        ier = LIB.mnt_regridedges_loadWeights(self.obj, fn, len(fn))
        if ier:
            error_handler(FILE, 'loadWeights', ier)

    def apply(self, srcdata, dstdata):
        """
        Apply the regridding weights to an edge field with unique edge Ids

        :param srcdata: contiguous arrays of source field data
        :param dstdata: contiguous arrays of destination field data (output)
        """
        LIB.mnt_regridedges_apply.argtypes = [POINTER(c_void_p),
                                              DOUBLE_ARRAY_PTR,
                                              DOUBLE_ARRAY_PTR]
        ier = LIB.mnt_regridedges_apply(self.obj, srcdata, dstdata)
        if ier:
            error_handler(FILE, 'apply', ier)
