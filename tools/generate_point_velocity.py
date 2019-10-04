import netCDF4
import numpy
import argparse
import sys
import time
import functools

parser = argparse.ArgumentParser(description='Generate point velocity')
parser.add_argument('-u', dest='velocityX', default='0.01*sin(pi*x/180.)',
                    help='Specify the contravariant velocity component u (deg/time) along longitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-v', dest='velocityY', default='0.01*cos(pi*y/180.)',
                    help='Specify the contravariant velocity component v (deg/time) along latitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid geometry/topology')
parser.add_argument('-N', dest='grid_var', default='', 
                    help='Specify the grid variable name in the netcdf file')
parser.add_argument('-d', dest='data_file', default='', 
                    help='Specify the netcdf file containing the edge integrated velocity data')
args = parser.parse_args()

if len(args.grid_file) == 0:
    print('ERROR: must provide grid file (-g)')
    sys.exit(1)

if len(args.data_file) == 0:
    print('ERROR: must provide data file (-d)')
    sys.exit(2)

if len(args.grid_var) == 0:
    print('ERROR: must provide a grid variable name (-N)')
    sys.exit(3)


nc = netCDF4.Dataset(args.grid_file, 'r')

# get the vertex coordinate names
xname, yname = nc[args.grid_var].node_coordinates.split(' ')

# read the coordinates
x = nc.variables[xname][:]
y = nc.variables[yname][:]

numPointDimName = nc.variables[xname].dimensions[0]

nc.close()

# write the velocity to disk
if args.data_file == args.grid_file:
    # append to grid file
    print('NOTE: point velocity field appended to grid file {}'.format(args.grid_file))
    nc = netCDF4.Dataset(args.data_file, 'r+')
else:
    # new file
    print('NOTE: point velocity field saved in {}'.format(args.data_file))
    nc = netCDF4.Dataset(args.data_file, 'w')
# write global attributes
nc.date = 'Created on {}'.format(time.asctime())
nc.command = functools.reduce(lambda x, y: x+' '+y, sys.argv)

numPoints = x.shape[0]
twoDimName = 'velocity_num_components'

# create dienmsions (unless there is already a dimension with the same name)
try:
    nc.createDimension(twoDimName, 2)
except:
    pass
try:
    nc.createDimension(numPointDimName, numPoints)
except:
    pass
velocity = nc.createVariable('velocity', 'f8', (numPointDimName, twoDimName))
velocity.location = 'node'

from numpy import sin, cos, pi
velocity[:, 0] = eval(args.velocityX)
velocity[:, 1] = eval(args.velocityY)
nc.close()

