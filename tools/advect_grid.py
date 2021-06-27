import numpy
import argparse
import sys
import time
from scipy.integrate import odeint
from ugrid_reader import UgridReader
import re

parser = argparse.ArgumentParser(description='Generate upstream grid')
parser.add_argument('-u', dest='velocity_x', default='sin(4*pi*(2*x-y)/180.)',
                    help='Specify the contravariant velocity component u (deg/time) along longitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-v', dest='velocity_y', default='cos(4*pi*(x-2.76*y)/180.)*cos(pi*y/180.)',
                    help='Specify the contravariant velocity component v (deg/time) along latitudes as a function of x (deg. east) and y (deg. north)')
parser.add_argument('-i', dest='input_grid_file', default='', 
                    help='Specify the netcdf file containing the grid geometry/topology and grid name as FILE_NAME:GRID_NAME')
parser.add_argument('-R', dest='regularization', action='store_false', 
                    help='Turn off grid regularization (recommended for uniform lat-lon)')
parser.add_argument('-t', dest='time', default=1.0, type=float,
                    help='Specify time step')
parser.add_argument('-n', dest='num_steps', default=1, type=int,
                    help='Specify number of steps')
parser.add_argument('-o', dest='output_vtk_files', default='', 
                    help='Specify the output VTK files for each time step')
args = parser.parse_args()

if len(args.input_grid_file) == 0:
    print('ERROR: must provide grid file (-i)')
    sys.exit(1)

if len(args.output_vtk_files) == 0:
    print('ERROR: must provide output VTK file (-o)')
    sys.exit(2)

ug = UgridReader(args.input_grid_file, regularization=args.regularization)
numCells = ug.getNumberOfCells()
print('number of cells = {}'.format(numCells))
numPoints = numCells * 4

vxy = numpy.zeros((numPoints*2,), numpy.float64)
xy =  numpy.zeros((numPoints*2,), numpy.float64)

# velocity field
from numpy import sin, cos, pi
def tendency(xy, t):
    x = xy[:numPoints]
    y = xy[numPoints:]
    vxy[:numPoints] = eval(args.velocity_x)
    vxy[numPoints:] = eval(args.velocity_y)
    return vxy

# advect grid forward
bname = re.sub(r'.vtk', '', args.output_vtk_files)
for i in range(args.num_steps):
    ug.saveToVtkFile(bname + "{:05d}".format(i) + '.vtk')
    x, y = ug.getLonLat()
    xy[:numPoints] = x.flat
    xy[numPoints:] = y.flat
    xyNew = odeint(tendency, xy, [0.0, args.time])
    ug.setLonLat(xyNew[1, :numPoints], xyNew[1, numPoints:])

ug.saveToVtkFile(bname + "{:05d}".format(i) + '.vtk')

