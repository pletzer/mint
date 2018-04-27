import numpy
import netCDF4
import argparse

from latlon import LatLon

parser = argparse.ArgumentParser()
parser.add_argument('-src', dest='src', type=str, default='latlon.nc', help='Source NetCDF file')
parser.add_argument('-vars', dest='vars', type=str, default='u,v', help='Name of u, v fields in source file')
parser.add_argument('-dst', dest='dst', type=str, default='', help='Destination NetCDF file')

args = parser.parse_args()

ll = LatLon()
ll.load(args.src)

