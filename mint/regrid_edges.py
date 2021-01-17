from ctypes import c_void_p, c_double, c_int, byref
import sys


# method signatures
MINTLIB.mnt_regridedges_new.argtypes = (byref(c_void_p),)
MINTLIB.mnt_regridedges_del.argtypes = (byref(c_void_p),)

class RegridEdges(object):


    def __init__(self):
        self.regrid_obj = c_void_p()
        ier = MINTLIB.mnt_regridedges_new(byref(self.regrid_obj))
        if ier: raise RuntimeError(f'ERROR: after calling {sys._getframe().f_code.co_filename} ier={ier}')


    def __del__(self):
        ier = MINTLIB.mnt_regridedges_del(byref(self.regrid_obj))
        if ier: raise RuntimeError(f'ERROR: after calling {sys._getframe().f_code.co_filename} ier={ier}')
        


    # def setSrcGridFlags(self, fixLonAcrossDateline, int averageLonAtPole):
    #     pass


    # def setDstGridFlags(self, fixLonAcrossDateline, int averageLonAtPole):
    #     pass


