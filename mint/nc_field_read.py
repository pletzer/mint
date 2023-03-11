from ctypes import (c_void_p, c_int, byref, POINTER, c_char_p,
                    c_size_t, create_string_buffer)
from . import MINTLIB
from . import error_handler
import numpy
import netCDF4
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

        self.nc = netCDF4.Dataset(fileName)
        self.var = self.nc.variables[varName]


    def getNumDims(self):
        """
        Get the number of dimensions.

        :return: number
        """
        return len(self.var.shape)


    def getDimName(self, iAxis):
        """
        Get the dimension name

        :param iAxis: axis index
        :return: name
        """
        return self.var.dimensions[iAxis]

    def getDim(self, iAxis):
        """
        Get the dimnension of an axis.

        :param iAxis: axis index
        :return: size
        """
        dimName = self.getDimName(iAxis)
        return len(self.nc.dimensions[dimName])

    def data(self, data=None):
        """
        Read the data.

        :param data: array to fill in (None if new array should be created)
        :return: return array
        """
        if data is None:
            return self.var[:]
        else:
            data[:] = self.var[:]
        return data


    def dataSlice(self, startInds0, counts, data=None):
        """
        Read a data slice.

        :param startInds0: start indices
        :param counts: counts along each dimension
        :param data: array to fill in (None if new array should be created)
        :return: return array
        """
        raise NotImplementedError('no yet implemented')
