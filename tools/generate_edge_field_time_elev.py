import netCDF4
import numpy
import argparse
import sys
import time
import functools

parser = argparse.ArgumentParser(description='Generate edge field')
parser.add_argument('-s', dest='streamFunction', default='sin(x*pi/180.)*cos(y*pi/180.)*z**2 * t',
                    help='Specify the stream function of x (deg. east), y (deg. north), z (level) and t (time)')
parser.add_argument('-nz', dest='nz', type=int, default=1,
                    help='Specify the number of vertical level cells')
parser.add_argument('-nt', dest='nt', type=int, default=1,
                    help='Specify the number of time intervals')
parser.add_argument('-n', dest='fieldName', default='edge_integrated_velocity',
                    help='Specify the name of the edge field')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='The netcdf file containing the grid geometry/topology and the name of the grid FILE_NAME$GRID_NAME')
parser.add_argument('-d', dest='data_file', default='', 
                    help='Specify the netcdf file containing the edge integrated velocity data')
args = parser.parse_args()

if len(args.grid_file) == 0:
    print('ERROR: must provide grid file (-g)')
    sys.exit(1)

if len(args.data_file) == 0:
    print('ERROR: must provide data file (-d)')
    sys.exit(2)

try:
    grid_file, grid_var = args.grid_file.split('$')
except:
    print('ERROR: could not extract grid name, specify -g FILE_NAME$GRID_NAME')
    sys.exit(3)


nc = netCDF4.Dataset(grid_file, 'r')

# get the vertex coordinate names
xname, yname = nc[grid_var].node_coordinates.split(' ')
edgeNodeConnName = nc[grid_var].edge_node_connectivity

# read the coordinates
x = nc.variables[xname][:]
y = nc.variables[yname][:]

# read the edge-node connectivity
edgeNodeConnectivity = nc.variables[edgeNodeConnName][:]

# subtract base index
start_index = getattr(nc.variables[edgeNodeConnName], 'start_index', 0)
edgeNodeConnectivity -= start_index

numEdgesDimName = nc.variables[edgeNodeConnName].dimensions[0]

nc.close()

# build elevations
zs = numpy.array(range(args.nz + 1), numpy.float64)

# build time values
ts = numpy.array(range(args.nt + 1), numpy.float64)

numPoints = x.shape[0]
numEdges = edgeNodeConnectivity.shape[0]

# compute the stream function on nodes
from numpy import sin, cos, pi, heaviside, power, exp
print(f'numEdges = {numEdges} x.shape = {x.shape}')
streamValues = numpy.zeros((args.nt + 1, args.nz + 1, numPoints), numpy.float64)
for l in range(args.nt + 1):
    t = ts[l]
    for k in range(args.nz + 1):
        z = zs[k]
        streamValues[l, k, :] = eval(args.streamFunction)

# compute the integrated velocity on edges
print(f'number of edges: {numEdges}')
integratedVelocity = numpy.zeros( (args.nt + 1, args.nz + 1, numEdges,), numpy.float64 )
i0, i1 = edgeNodeConnectivity[:, 0], edgeNodeConnectivity[:, 1]

for l in range(args.nt + 1):
    for k in range(args.nz + 1):
        integratedVelocity[l, k, :] = streamValues[l, k, i1] - streamValues[l, k, i0]

# write the velocity to disk
if args.data_file == grid_file:
    # append to grid file
    print(f'note: edge field "{args.fieldName}" will be appended to ugrid {grid_file}:{grid_var}')
    nc = netCDF4.Dataset(args.data_file, 'r+')
else:
    # new file
    print('note: saving edge field "{}" in {}'.format(args.fieldName, args.data_file))
    nc = netCDF4.Dataset(args.data_file, 'w')
# write global attributes
nc.date = 'Created on {}'.format(time.asctime())
nc.command = functools.reduce(lambda x, y: x+' '+y, sys.argv)
try:
    nc.createDimension(numEdgesDimName, numEdges)
except:
    # might already have a dimension called numEdgesDimName in the file
    pass

numElevDimName = 'nelevs'
numTimeDimName = 'ntimes'
nc.createDimension(numElevDimName, args.nz + 1)
nc.createDimension(numTimeDimName, args.nt + 1)
vel = nc.createVariable(args.fieldName, 'f8', (numTimeDimName, numElevDimName, numEdgesDimName,))
vel.location = 'edge'
vel.mesh = grid_var
vel[:] = integratedVelocity
nc.close()

