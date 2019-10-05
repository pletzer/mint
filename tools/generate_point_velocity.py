import netCDF4
import numpy
import argparse
import sys
import time
import functools
from scipy.integrate import odeint


parser = argparse.ArgumentParser(description='Generate point velocity')
parser.add_argument('-u', dest='velocityX', default='0.01*sin(pi*x/180.)',
                    help='Specify the contravariant velocity component u (deg/time) along longitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-v', dest='velocityY', default='0.01*cos(pi*y/180.)',
                    help='Specify the contravariant velocity component v (deg/time) along latitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid geometry/topology')
parser.add_argument('-N', dest='grid_var', default='', 
                    help='Specify the grid variable name in the netcdf file')
parser.add_argument('-tf', dest='finalTime', default=1.0, type=float,
                    help='Specify final time for integrating trajectories upstream')
parser.add_argument('-d', dest='upstream_data_file', default='', 
                    help='Specify the netcdf file containing the upstream coordinates and the velocity')
args = parser.parse_args()

if len(args.grid_file) == 0:
    print('ERROR: must provide grid file (-g)')
    sys.exit(1)

if len(args.upstream_data_file) == 0:
    print('ERROR: must provide upstream data file (-d)')
    sys.exit(2)

if args.upstream_data_file == args.grid_file:
    print('ERROR: upstream data file name must be different from grid file name')
    sys.exit(3)

if len(args.grid_var) == 0:
    print('ERROR: must provide a grid variable name (-N)')
    sys.exit(4)


def copyAttributes(fromNcVar, toNcVar):
    """
    Copy the attributes from fromNcVar to toNcVar
    """
    for attrName in fromNcVar.ncattrs():
        attrVal = getattr(fromNcVar, attrName)
        setattr(toNcVar, attrName, attrVal)


ncGrid = netCDF4.Dataset(args.grid_file, 'r')

# get the vertex coordinate names
xName, yName = ncGrid[args.grid_var].node_coordinates.split(' ')

# get grid dimensions
numPoints = ncGrid[xName].shape[0]
numEdges = ncGrid[ncGrid[args.grid_var].edge_node_connectivity].shape[0]
numFaces = ncGrid[ncGrid[args.grid_var].face_node_connectivity].shape[0]

# read the coordinates (initial conditions)
x = ncGrid.variables[xName][:]
y = ncGrid.variables[yName][:]

# open the output file
ncUp = netCDF4.Dataset(args.upstream_data_file, 'w')
# write global attributes
ncUp.date = 'Created on {}'.format(time.asctime())
ncUp.command = functools.reduce(lambda x, y: x+' '+y, sys.argv)

# create dimensions and variables
twoDimName = 'two'
fourDimName = 'four'
numPointsDimName = 'num_points'
numEdgesDimName = 'num_edges'
numFacesDimName = 'num_faces'
ncUp.createDimension(twoDimName, 2)
ncUp.createDimension(fourDimName, 4)
ncUp.createDimension(numPointsDimName, numPoints)
ncUp.createDimension(numEdgesDimName, numEdges)
ncUp.createDimension(numFacesDimName, numFaces)

# copy the topology over
gridVarUpName = args.grid_var + '_upstream'
gridVarUp = ncUp.createVariable(gridVarUpName, 'i4')
# save upstream grid
for attrName in ncGrid[args.grid_var].ncattrs():
    attrVal = getattr(ncGrid[args.grid_var], attrName)
    setattr(gridVarUp, attrName, attrVal)
# new coordinates
xNameUp = xName + '_upstream'
yNameUp = yName + '_upstream'
gridVarUp.node_coordinates = '{} {}'.format(xNameUp, yNameUp)

faceNodeUp = ncUp.createVariable(gridVarUp.face_node_connectivity, 'i4', 
                               (numFacesDimName, fourDimName))
copyAttributes(ncGrid[ncGrid[args.grid_var].face_node_connectivity], faceNodeUp)
faceNodeUp[:] = ncGrid[ncGrid[args.grid_var].face_node_connectivity][:]

faceEdgeUp = ncUp.createVariable(gridVarUp.face_edge_connectivity, 'i4',
                               (numFacesDimName, fourDimName))
copyAttributes(ncGrid[ncGrid[args.grid_var].face_edge_connectivity], faceEdgeUp)
faceEdgeUp[:] = ncGrid[ncGrid[args.grid_var].face_edge_connectivity][:]


edgeNodeUp = ncUp.createVariable(gridVarUp.edge_node_connectivity, 'i4',
                               (numEdgesDimName, twoDimName))
copyAttributes(ncGrid[ncGrid[args.grid_var].edge_node_connectivity], edgeNodeUp)
edgeNodeUp[:] = ncGrid[ncGrid[args.grid_var].edge_node_connectivity][:]


velocity = ncUp.createVariable('velocity', 'f8', (numPointsDimName, twoDimName))
velocity.location = 'node'

# velocity at nodes
from numpy import sin, cos, pi
velocity[:, 0] = eval(args.velocityX)
velocity[:, 1] = eval(args.velocityY)

# integrate the nodal positions backwards in time
vxy = numpy.zeros((numPoints*2,), numpy.float64)
def tendency(xy, t):
    x = xy[:numPoints]
    y = xy[numPoints:]
    vxy[:numPoints] = eval(args.velocityX)
    vxy[numPoints:] = eval(args.velocityY)
    return vxy

xy = numpy.zeros((numPoints*2,), numpy.float64)
xy[:numPoints] = x
xy[numPoints:] = y
xyUpstream = odeint(tendency, xy, [0.0, -args.finalTime])
xUpstream = xyUpstream[1, :numPoints]
yUpstream = xyUpstream[1, numPoints:]


# save the new coordinates
xVarUp = ncUp.createVariable(xNameUp, 'f8', (numPointsDimName,))
copyAttributes(ncGrid[xName], xVarUp)

yVarUp = ncUp.createVariable(yNameUp, 'f8', (numPointsDimName,))
copyAttributes(ncGrid[yName], yVarUp)

# write
xVarUp[:] = xUpstream
yVarUp[:] = yUpstream

ncGrid.close()
ncUp.close()

