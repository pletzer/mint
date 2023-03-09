from ctypes import (c_void_p, c_int, byref, POINTER, c_char_p,
                    c_size_t, create_string_buffer)
from . import MINTLIB
from . import error_handler
import numpy
import logging


FILE = 'nc_field_read.py'
DOUBLE_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.float64)
SIZET_ARRAY_PTR = numpy.ctypeslib.ndpointer(dtype=numpy.uintp)
MAX_DIM_NAME = 1024

MINTLIB.nc_open.argtypes = [c_char_p, c_int, POINTER(c_int)]
MINTLIB.nc_close.argtypes = [c_int]


class NcFieldRead(object):
    """
    A class to read edge fields stored in NetCDF files.
    """

    def __init__(self, fileName, varName):
        """
        Constructor.

        :param fileName: NetCDF file name
        :param varName: variable name
        """

        self.ptr = c_void_p()
        self.obj = byref(self.ptr)

        # open the netcdf file and hang on to the handle
        self.ncid = c_int(-1)
        ier = MINTLIB.nc_open(str(fileName).encode('utf-8'), 0, byref(self.ncid))
        if ier != 0:
            # file does not exist?
            logging.error(f'file {fileName} could not be opened')
            return

        # get the variable Id
        self.varid = c_int(-1)
        ier = MINTLIB.nc_inq_varid(self.ncid.value, varName.encode('utf-8'), byref(self.varid))
        if ier != 0:
            # variable does not exist?
            logging.error(f'could not get the var Id of {varName} in file {fileName}')
            ier = MINTLIB.nc_close(self.ncid.value)
            return

        MINTLIB.mnt_ncfieldread_new.argtypes = [POINTER(c_void_p),
                                                c_int, c_int]
        ier = MINTLIB.mnt_ncfieldread_new(self.obj, self.ncid, self.varid)
        if ier:
            error_handler(FILE, '__init__', ier)

    def __del__(self):
        """
        Destructor.
        """
        MINTLIB.mnt_ncfieldread_del.argtypes = [POINTER(c_void_p)]
        ier = MINTLIB.mnt_ncfieldread_del(self.obj)
        if ier:
            error_handler(FILE, '__del__', ier)

        # close the file
        if self.ncid.value >= 0:
            MINTLIB.nc_close(self.ncid.value)


    def getNumDims(self):
        """
        Get the number of dimensions.

        :return: number
        """
        MINTLIB.mnt_ncfieldread_getNumDims.argtypes = [POINTER(c_void_p),
                                                       POINTER(c_int)]
        size = c_int()
        ier = MINTLIB.mnt_ncfieldread_getNumDims(self.obj, byref(size))
        if ier:
            error_handler(FILE, 'getNumDims', ier)
        return size.value

    def getDimName(self, iAxis):
        """
        Get the dimension name

        :param iAxis: axis index
        :return: name
        """
        MINTLIB.mnt_ncfieldread_getDimName.argtypes = [POINTER(c_void_p),
                                                       c_int, c_char_p,
                                                       c_int]
        dimName = b" "*MAX_DIM_NAME
        dimNameSize = len(dimName)
        ier = MINTLIB.mnt_ncfieldread_getDimName(self.obj,
                                                 iAxis,
                                                 dimName,
                                                 dimNameSize)
        if ier:
            error_handler(FILE, 'getDimName', ier)
        return dimName.rstrip(b'\x00')

    def getDim(self, iAxis):
        """
        Get the dimnension of an axis.

        :param iAxis: axis index
        :return: size
        """
        MINTLIB.mnt_ncfieldread_getDim.argtypes = [POINTER(c_void_p),
                                                   c_int, POINTER(c_size_t)]
        size = c_size_t()
        ier = MINTLIB.mnt_ncfieldread_getDim(self.obj, iAxis, byref(size))
        if ier:
            error_handler(FILE, 'getDim', ier)
        return size.value

    def data(self, data=None):
        """
        Read the data.

        :param data: array to fill in (None if new array should be created)
        :return: return array
        """
        MINTLIB.mnt_ncfieldread_data.argtypes = [POINTER(c_void_p),
                                                 DOUBLE_ARRAY_PTR]
        if data is None:
            dims = [self.getDim(iAxis) for iAxis in range(self.getNumDims())]
            data = numpy.empty(dims, numpy.float64)
        ier = MINTLIB.mnt_ncfieldread_data(self.obj, data)
        if ier:
            error_handler(FILE, 'data', ier)
        return data

    def dataSlice(self, startInds0, counts, data=None):
        """
        Read a data slice.

        :param startInds0: start indices
        :param counts: counts along each dimension
        :param data: array to fill in (None if new array should be created)
        :return: return array
        """
        MINTLIB.mnt_ncfieldread_dataSlice.argtypes = [POINTER(c_void_p),
                                                      SIZET_ARRAY_PTR,
                                                      SIZET_ARRAY_PTR,
                                                      DOUBLE_ARRAY_PTR]
        if data is None:
            data = numpy.empty(counts, numpy.float64)
        ier = MINTLIB.mnt_ncfieldread_dataSlice(self.obj, startInds0,
                                                counts, data)
        if ier:
            error_handler(FILE, 'dataSlice', ier)
        return data
