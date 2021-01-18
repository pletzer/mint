from ctypes import c_void_p, c_double, c_int, byref, POINTER, c_char_p
from . import LIB

def error_handler(filename, methodname, ier):
    raise RuntimeError(f'ERROR ier={ier} after calling {methodname} in {filename}!')

class RegridEdges(object):


    def __init__(self):
        """
        Regrid edge field constructor.
        """

        self.regrid_obj = byref(c_void_p())

        ier = LIB.mnt_regridedges_new(self.regrid_obj)
        if ier: 
            error_handler('regrid_edges.py', '__init__', ier)


    def __del__(self):
        """
        Regrid edge field destructor.
        """
        ier = LIB.mnt_regridedges_del(self.regrid_obj)
        if ier: 
            error_handler('regrid_edges.py', '__del__', ier)


    def setSrcGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the source grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length should be added/subtracted
                                     in order to make eeach cell as compact as poassible
        :param averageLonAtPole: set to 1 if the longitudes at the poles should the average of 
                                 the cell's longitudes

        note:: a lon-lat grid requires 0, 0 and a cibed sphere grid requires 1, 1
        """
        LIB.mnt_regridedges_setSrcGridFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = LIB.mnt_regridedges_setSrcGridFlags(self.regrid_obj, fixLonAcrossDateline, averageLonAtPole)
        if ier:
            error_handler('regrid_edges.py', 'setSrcGridFlags', ier)


    def setDstGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        """
        Set the destrination grid flags.

        :param fixLonAcrossDateline: set to 1 if a periodicity length should be added/subtracted
                                     in order to make eeach cell as compact as poassible
        :param averageLonAtPole: set to 1 if the longitudes at the poles should the average of 
                                 the cell's longitudes

        note:: a lon-lat grid requires 0, 0 and a cibed sphere grid requires 1, 1
        """
        LIB.mnt_regridedges_setDstGridFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = LIB.mnt_regridedges_setDstGridFlags(self.regrid_obj, fixLonAcrossDateline, averageLonAtPole)
        if ier:
            error_handler('regrid_edges.py', 'setDstGridFlags', ier)


    def loadSrcGrid(self, filename):
        """
        Load a source grid from a 2D UGRID file.

        :param filename: string in the format filename:meshname
        """
        LIB.mnt_regridedges_loadSrcGrid.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadSrcGrid(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadSrcGrid', ier)


    def loadDstGrid(self, filename):
        """
        Load a destination grid from a 2D UGRID file.

        :param filename: string in the format filename:meshname
        """
        LIB.mnt_regridedges_loadDstGrid.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadDstGrid(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadDstGrid', ier)


    def build(self, numCellsPerBucket, periodX, debug):
        """
        Build the regridder and compute the regridding weights

        :param numCellsPerBucket: average number of cells per bucket (affects performance only)
        :param periodX: periodicity length (set to 0 if non-periodic)
        :param debug: 0=no debug info, 1=print debug info, 2=save bad edges in VTK file
        """
        LIB.mnt_regridedges_build.argtypes = [POINTER(c_void_p), c_int, c_double, c_int]
        ier = LIB.mnt_regridedges_build(self.regrid_obj, numCellsPerBucket, periodX, debug)
        if ier:
            error_handler('regrid_edges.py', 'build', ier)


    def dumpWeights(self, filename):
        """
        Dump the weights to a file
        :param filename: file name
        """
        LIB.mnt_regridedges_dumpWeights.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_dumpWeights(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'dumpWeights', ier)


    def loadWeights(self, filename):
        """
        Load the weights from a file
        :param filename: file name
        """
        LIB.mnt_regridedges_loadWeights.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadWeights(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadWeights', ier)


    def apply(self, srcdata, dstdata):
        """
        Apply the regridding weights to an edge field with unique edge Ids

        :param srcdata: contiguous arrays of source field data
        :param dstdata: contiguous arrays of destination field data (will be filled in)
        """
        LIB.mnt_regridedges_apply.argtypes = [POINTER(c_void_p), 
                                              numpy.ctypeslib.ndpointer(dtype=numpy.float64), 
                                              numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        ier = LIB.mnt_regridedges_apply(srcdata, dstdata)
        if ier:
            error_handler('regrid_edges.py', 'apply', ier)


