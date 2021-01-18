from ctypes import c_void_p, c_double, c_int, byref, POINTER, c_char_p
from . import LIB

def error_handler(filename, methodname, ier):
    raise RuntimeError(f'ERROR ier={ier} after calling {methodname} in {filename}!')

class RegridEdges(object):


    def __init__(self):

        self.regrid_obj = byref(c_void_p())

        ier = LIB.mnt_regridedges_new(self.regrid_obj)
        if ier: 
            error_handler('regrid_edges.py', '__init__', ier)


    def __del__(self):
        ier = LIB.mnt_regridedges_del(self.regrid_obj)
        if ier: 
            error_handler('regrid_edges.py', '__del__', ier)


    def setSrcGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        LIB.mnt_regridedges_setSrcGridFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = LIB.mnt_regridedges_setSrcGridFlags(self.regrid_obj, fixLonAcrossDateline, averageLonAtPole)
        if ier:
            error_handler('regrid_edges.py', 'setSrcGridFlags', ier)


    def setDstGridFlags(self, fixLonAcrossDateline, averageLonAtPole):
        LIB.mnt_regridedges_setDstGridFlags.argtypes = [POINTER(c_void_p), c_int, c_int]
        ier = LIB.mnt_regridedges_setDstGridFlags(self.regrid_obj, fixLonAcrossDateline, averageLonAtPole)
        if ier:
            error_handler('regrid_edges.py', 'setDstGridFlags', ier)


    def loadSrcGrid(self, filename):
        LIB.mnt_regridedges_loadSrcGrid.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadSrcGrid(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadSrcGrid', ier)


    def loadDstGrid(self, filename):
        LIB.mnt_regridedges_loadDstGrid.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadDstGrid(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadDstGrid', ier)


    def build(self, numCellsPerBucket, periodX, debug):
        LIB.mnt_regridedges_build.argtypes = [POINTER(c_void_p), c_int, c_double, c_int]
        ier = LIB.mnt_regridedges_build(self.regrid_obj, numCellsPerBucket, periodX, debug)
        if ier:
            error_handler('regrid_edges.py', 'build', ier)


    def dumpWeights(self, filename):
        LIB.mnt_regridedges_dumpWeights.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_dumpWeights(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'dumpWeights', ier)


    def loadWeights(self, filename):
        LIB.mnt_regridedges_loadWeights.argtypes = [POINTER(c_void_p), c_char_p, c_int]
        ier = LIB.mnt_regridedges_loadWeights(self.regrid_obj, filename, len(filename))
        if ier:
            error_handler('regrid_edges.py', 'loadWeights', ier)


    def apply(self, srcdata, dstdata):
        LIB.mnt_regridedges_apply.argtypes = [POINTER(c_void_p), 
                                              numpy.ctypeslib.ndpointer(dtype=numpy.float64), 
                                              numpy.ctypeslib.ndpointer(dtype=numpy.float64)]
        ier = LIB.mnt_regridedges_apply(srcdata, dstdata)
        if ier:
            error_handler('regrid_edges.py', 'apply', ier)


