from ctypes import (c_void_p, c_int, byref, POINTER, c_char_p,
                    c_size_t, c_longlong, c_double)
from . import MINTLIB, NUM_VERTS_PER_QUAD, NUM_VERTS_PER_EDGE
from . import error_handler
import numpy


FILE = 'multi_array_iter.py'
SIZE_T_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.uintp)


class MultiArrayIter(object):
    """
    A class to iterate data across time, elevation and other non-spatial dimensions
    """

    def __init__(self, dims):
        """
        Constructor.

        :param dims array of data dimensions
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        # store the multi-indices
        self.indices = numpy.empty(len(dims), numpy.uintp)

        MINTLIB.mnt_multiarrayiter_new.argtypes = [POINTER(c_void_p), c_int, SIZE_T_ARRAY_PTR]
        ier = MINTLIB.mnt_multiarrayiter_new(self.obj, len(dims), numpy.array(dims, numpy.uint64))
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_multiarrayiter_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_multiarrayiter_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

    def begin(self):
        """
        Set the iterator to the beginning.
        """
        MINTLIB.mnt_multiarrayiter_begin.argtypes = [POINTER(c_void_p),]
        ier = MINTLIB.mnt_multiarrayiter_begin(self.obj)
        if ier:
            error_handler(FILE, 'begin', ier)

    def next(self):
        """
        Increment the iterator.

        :note: no test will be performed if the iterator steps out of bounds.
        """
        MINTLIB.mnt_multiarrayiter_next.argtypes = [POINTER(c_void_p),]
        ier = MINTLIB.mnt_multiarrayiter_next(self.obj)
        if ier:
            error_handler(FILE, 'next', ier)


    def getNumIters(self):
        """
        Get the number of iterations.

        :returns number
        """
        MINTLIB.mnt_multiarrayiter_getNumIters.argtypes = [POINTER(c_void_p), POINTER(c_size_t)]
        n = c_size_t()
        ier = MINTLIB.mnt_multiarrayiter_getNumIters(self.obj, byref(n))
        if ier:
            error_handler(FILE, 'getNumIters', ier)
        return n.value

    def getIndices(self):
        """
        Get the indices at this iteration.

        :returns array of indices
        """
        MINTLIB.mnt_multiarrayiter_getIndices.argtypes = [POINTER(c_void_p), SIZE_T_ARRAY_PTR]
        ier = MINTLIB.mnt_multiarrayiter_getIndices(self.obj, self.indices)
        if ier:
            error_handler(FILE, 'getIndices', ier)
        return self.indices

