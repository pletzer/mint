from ctypes import (c_char_p, c_size_t)
from . import MINTLIB
import logging

def error_handler(filename, methodname, ier, detailedmsg=''):
    """
    Error handler
    @param filename file name
    @param methodname method or function name
    @param ier error code
    """
    msg = f'ier={ier} after calling {methodname} in {filename}!\n'
    msg += detailedmsg
    logging.error(msg)
    raise RuntimeError(msg)

def warning_handler(filename, methodname, ier, detailedmsg=''):
    """
    Warning handler
    @param filename file name
    @param methodname method or function name
    @param ier error code
    """
    msg = f'ier={ier} after calling {methodname} in {filename}!\n'
    msg += detailedmsg
    logging.warning(msg)

def printLogMessages():
    """
    Print the log messages from the MINT library
    """
    MINTLIB.mnt_printLogMessages.restype = None
    MINTLIB.mnt_printLogMessages()

def writeLogMessages(filename):
    """
    Write the log messages from the MINT library to file
    @param filename output file name
    """
    MINTLIB.mnt_writeLogMessages.argtypes = [c_char_p, c_size_t]
    MINTLIB.mnt_writeLogMessages.restype = None
    fn = filename.encode('utf-8')
    MINTLIB.mnt_writeLogMessages(fn, len(fn))



