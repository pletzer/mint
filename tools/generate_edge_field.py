import netCDF4
import numpy
import argparse
import sys
import time
import functools

parser = argparse.ArgumentParser(description='Generate edge field')
parser.add_argument('-s', dest='streamFunction', default='sin(x*pi/180.)*cos(y*pi/180.)',
                    help='Specify the stream function of x (deg. east) and y (deg. north)')
parser.add_argument('-n', dest='fieldName', default='edge_integrated_velocity',
                    help='Specify the name of the edge field')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid geometry/topology')
parser.add_argument('-d', dest='data_file', default='', 
                    help='Specify the netcdf file containing the edge integrated velocity data')
args = parser.parse_args()

if len(args.grid_file) == 0:
    print('ERROR: must provide grid file (-g)')
    sys.exit(1)

if len(args.data_file) == 0:
    print('ERROR: must provide data file (-d)')
    sys.exit(2)

nc = netCDF4.Dataset(args.grid_file, 'r')

# read the coordinates
x = nc.variables['physics_node_x'][:]
y = nc.variables['physics_node_y'][:]

# read the edge-node connectivity
edgeNodeConnectivity = nc.variables['physics_edge_nodes'][:]

# subtract base index
start_index = getattr(nc.variables['physics_edge_nodes'], 'start_index', 0)
edgeNodeConnectivity -= start_index

numEdgesDimName = nc.variables['physics_edge_nodes'].dimensions[0]

nc.close()

# compute the stream function on nodes
from numpy import sin, cos, pi 
streamValues = eval(args.streamFunction)

# compute the integrated velocity on edges
numEdges = edgeNodeConnectivity.shape[0]
integratedVelocity = numpy.zeros( (numEdges,), numpy.float64 )
i0, i1 = edgeNodeConnectivity[:, 0], edgeNodeConnectivity[:, 1]
integratedVelocity = streamValues[i1] - streamValues[i0]

# write the velocity to disk
if args.data_file == args.grid_file:
    # append to grid file
    print('NOTE: edge field will be appended to grid file {}'.format(args.grid_file))
    nc = netCDF4.Dataset(args.data_file, 'r+')
else:
    # new file
    print('NOTE: saving edge field in {}'.format(args.data_file))
    nc = netCDF4.Dataset(args.data_file, 'w')
# write global attributes
nc.date = 'Created on {}'.format(time.asctime())
nc.command = functools.reduce(lambda x, y: x+' '+y, sys.argv)
try:
    nc.createDimension(numEdgesDimName, numEdges)
except:
    # might already have a dimension called numEdgesDimName in the file
    pass
vel = nc.createVariable(args.fieldName, 'f8', (numEdgesDimName,))
vel.location = 'edge'
vel[:] = integratedVelocity
nc.close()

