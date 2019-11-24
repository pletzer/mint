import netCDF4
import numpy
import argparse
import sys
import time
import functools
from scipy.integrate import odeint


parser = argparse.ArgumentParser(description='Generate upstream grid')
parser.add_argument('-u', dest='velocity_x', default='sin(pi*x/180.)',
                    help='Specify the contravariant velocity component u (deg/time) along longitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-v', dest='velocity_y', default='cos(pi*y/180.)',
                    help='Specify the contravariant velocity component v (deg/time) along latitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid geometry/topology and grid name as FILE_NAME:GRID_NAME')
parser.add_argument('-t', dest='time', default=1.0, type=float,
                    help='Specify final time for integrating trajectories upstream')
parser.add_argument('-o', dest='output_file', default='', 
                    help='Specify the output netcdf file containing the upstream coordinates and the velocity')
parser.add_argument('-P', dest='periodicityLength', default=360., 
                    help='Periodicity length in x (set to zero if non-periodic)')
args = parser.parse_args()

if len(args.grid_file) == 0:
    print('ERROR: must provide grid file (-g)')
    sys.exit(1)

if len(args.output_file) == 0:
    print('ERROR: must provide output file (-o)')
    sys.exit(2)

if args.output_file == args.grid_file:
    print('ERROR: output file name must be different from grid file name')
    sys.exit(3)

# file_name:grid_name
try:
    grid_file, grid_var = args.grid_file.split(':')
except:
    print('ERROR: argument passed to -g should be in format GRID_FILE:GRID_NAME')
    sys.exit(4)


def copyAttributes(fromNcVar, toNcVar):
    """
    Copy the attributes from fromNcVar to toNcVar
    """
    for attrName in fromNcVar.ncattrs():
        attrVal = getattr(fromNcVar, attrName)
        setattr(toNcVar, attrName, attrVal)


ncGrid = netCDF4.Dataset(grid_file, 'r')

try:
    ncGridVar = ncGrid[grid_var]
except:
    print('ERROR: could not extract grid name {} from file {}'.format(grid_var, grid_file))
    print('possible grid names are:')
    count = 0
    for v in ncGrid.variables:
        var = ncGrid[v]
        if hasattr(var, 'cf_role') and var.cf_role == 'mesh_topology':
            print(v)
            count += 1
    print('found {} variables with attribute "cf_role" == "mesh_topology"'.format(count))
    sys.exit(4)

# get the vertex coordinate names
xName, yName = ncGridVar.node_coordinates.split(' ')

# get grid dimensions
numPoints = ncGrid[xName].shape[0]
numEdges = ncGrid[ncGridVar.edge_node_connectivity].shape[0]
numFaces = ncGrid[ncGridVar.face_node_connectivity].shape[0]

# read the coordinates (initial conditions)
xInitial = ncGrid.variables[xName][:]
yInitial = ncGrid.variables[yName][:]

# open the output file
ncUp = netCDF4.Dataset(args.output_file, 'w')
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
gridVarUpName = grid_var
gridVarUp = ncUp.createVariable(gridVarUpName, 'i4')
# save upstream grid
for attrName in ncGridVar.ncattrs():
    attrVal = getattr(ncGridVar, attrName)
    setattr(gridVarUp, attrName, attrVal)
# new coordinates
xNameUp = xName
yNameUp = yName
gridVarUp.node_coordinates = '{} {}'.format(xNameUp, yNameUp)

faceNodeUp = ncUp.createVariable(gridVarUp.face_node_connectivity, 'i4', 
                               (numFacesDimName, fourDimName))
copyAttributes(ncGrid[ncGridVar.face_node_connectivity], faceNodeUp)
faceNodeUp[:] = ncGrid[ncGridVar.face_node_connectivity][:]

faceEdgeUp = ncUp.createVariable(gridVarUp.face_edge_connectivity, 'i4',
                               (numFacesDimName, fourDimName))
copyAttributes(ncGrid[ncGridVar.face_edge_connectivity], faceEdgeUp)
faceEdgeUp[:] = ncGrid[ncGridVar.face_edge_connectivity][:]


edgeNodeUp = ncUp.createVariable(gridVarUp.edge_node_connectivity, 'i4',
                               (numEdgesDimName, twoDimName))
copyAttributes(ncGrid[ncGridVar.edge_node_connectivity], edgeNodeUp)
edgeNodeUp[:] = ncGrid[ncGridVar.edge_node_connectivity][:]


velocity = ncUp.createVariable('velocity', 'f8', (numPointsDimName, twoDimName))
velocity.location = 'node'
velocity.mesh = grid_var

# velocity at nodes
from numpy import sin, cos, pi
x = xInitial
y = yInitial
velocity[:, 0] = eval(args.velocity_x)
velocity[:, 1] = eval(args.velocity_y)

# integrate the nodal positions backwards in time
vxy = numpy.zeros((numPoints*2,), numpy.float64)
def tendency(xy, t):
    x = xy[:numPoints]
    y = xy[numPoints:]
    vxy[:numPoints] = eval(args.velocity_x)
    vxy[numPoints:] = eval(args.velocity_y)
    return vxy

# integrate the trajectories upstream 
xy = numpy.zeros((numPoints*2,), numpy.float64)
xy[:numPoints] = xInitial
xy[numPoints:] = yInitial
xMin = xInitial.min()
xMax = xInitial.max()
xyUpstream = odeint(tendency, xy, [0.0, -args.time])

# we're only interested in the final positions
xUpstream = xyUpstream[1, :numPoints]
yUpstream = xyUpstream[1, numPoints:]

# apply periodicity
#xUpstream += (xUpstream < 0.) * args.periodicityLength
#xUpstream -= (xUpstream > args.periodicityLength) * args.periodicityLength

# make sure the latitudes are within [-90, 90]
numpy.clip(yUpstream, -90., 90., out=yUpstream)

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

