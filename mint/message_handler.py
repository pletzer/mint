import logging

def error_handler(filename, methodname, ier):
    msg = f'ier={ier} after calling {methodname} in {filename}!'
    logging.error(msg)
    raise RuntimeError(msg)

def warning_handler(filename, methodname, ier, detailedmsg=''):
    msg = f'ier={ier} after calling {methodname} in {filename}!\n'
    msg += detailedmsg
    logging.warning(msg)
