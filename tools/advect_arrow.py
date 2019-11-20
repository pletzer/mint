import vtk
import numpy
import sys
import argparse
from scipy.integrate import odeint

from math import sin, cos, pi
from numpy import linspace

"""
Generate grid and edge data on uniform grid and save result in UGRID file
"""

parser = argparse.ArgumentParser(description='Advect arrow')
parser.add_argument('-o', dest='grid_file', default='grid.vtk', 
                    help='Specify the output VTK file name containing the grid, stream function and velocity data')
parser.add_argument('-nx', default=1, type=int, 
                    help='Number of longitude cells')
parser.add_argument('-ny', default=1, type=int, 
                    help='Number of latitude cells')
parser.add_argument('-xline', default="linspace(160,160,101)", type=str, dest="initXPoints",
                    help='Initial x line coordinates')
parser.add_argument('-yline', default="linspace(10,80,101)", type=str, dest="initYPoints",
                    help='Initial y line coordinates')
parser.add_argument('-s', type=str, dest='stream_funct', default='100*(sin((pi*(x - 2*y))/360.) + sin((pi*(x + 2*y))/360.))/2.', 
                   help='Stream function of x (longitude in deg) and y (latitude in deg) used for setting the edge integrals')
parser.add_argument('-u', type=str, dest='u_funct', default='100*((pi*cos((pi*(x - 2*y))/360.))/180. - (pi*cos((pi*(x + 2*y))/360.))/180.)/2.', 
                   help='u function of x (longitude in deg) and y (latitude in deg)')
parser.add_argument('-v', type=str, dest='v_funct', default='100*((pi*cos((pi*(x - 2*y))/360.))/360. + (pi*cos((pi*(x + 2*y))/360.))/360.)/2.', 
                   help='v function of x (longitude in deg) and y (latitude in deg)')
parser.add_argument('-t', default=1.0, type=float, dest="timeStep",
                    help='Time step')
parser.add_argument('-nt', default=10, type=int, dest="numSteps",
                    help='Number of time steps')
args = parser.parse_args()

# check
if not args.grid_file:
    print("ERROR: must specify grid file (-o)")
    sys.exit(1)


nx, ny = args.nx, args.ny
nnodes = (nx + 1) * (ny + 1)
nedges = nx * (ny + 1) + (nx + 1) * ny
nfaces = nx * ny
npoints = nfaces * 4

dx = 360.0/float(nx)
dy = 180.0/float(ny)

pdata = vtk.vtkDoubleArray()
points = vtk.vtkPoints()
grid = vtk.vtkUnstructuredGrid()
sdata = vtk.vtkDoubleArray()
vdata = vtk.vtkDoubleArray()
writer = vtk.vtkUnstructuredGridWriter()

pdata.SetNumberOfComponents(3)
pdata.SetNumberOfTuples(npoints)

sdata.SetNumberOfComponents(1)
sdata.SetNumberOfTuples(npoints)
sdata.SetName('stream_function')

vdata.SetNumberOfComponents(3)
vdata.SetNumberOfTuples(nfaces)
vdata.SetName('velocity')


iface = 0
for j in range(ny):
    y0 = -90.0 + j*dy
    y1 = y0 + dy
    for i in range(nx):
        x0 = 0.0 + i*dx
        x1 = x0 + dx
        pdata.SetTuple(iface*4 + 0, (x0, y0, 0.))
        pdata.SetTuple(iface*4 + 1, (x1, y0, 0.))
        pdata.SetTuple(iface*4 + 2, (x1, y1, 0.))
        pdata.SetTuple(iface*4 + 3, (x0, y1, 0.))
        x, y = x0, y0
        s0 = eval(args.stream_funct)
        x, y = x1, y0
        s1 = eval(args.stream_funct)
        x, y = x1, y1
        s2 = eval(args.stream_funct)
        x, y = x0, y1
        s3 = eval(args.stream_funct)
        sdata.SetTuple(iface*4 + 0, (s0,))
        sdata.SetTuple(iface*4 + 1, (s1,))
        sdata.SetTuple(iface*4 + 2, (s2,))
        sdata.SetTuple(iface*4 + 3, (s3,))
        x = 0.5*(x0 + x1)
        y = 0.5*(y0 + y1)
        u = eval(args.u_funct)
        v = eval(args.v_funct)
        vdata.SetTuple(iface, (u, v, 0.0))
        iface += 1

points.SetData(pdata)
grid.SetPoints(points)
grid.Allocate(nfaces, 1)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)
for iface in range(nfaces):
    ptIds.SetId(0, iface*4 + 0)
    ptIds.SetId(1, iface*4 + 1)
    ptIds.SetId(2, iface*4 + 2)
    ptIds.SetId(3, iface*4 + 3)
    grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

grid.GetPointData().AddArray(sdata)
grid.GetCellData().AddArray(vdata)

# write to file
print(f'writing grid to {args.grid_file}')
writer.SetFileName(args.grid_file)
writer.SetInputData(grid)
writer.Update()

initXPoints = eval(args.initXPoints)
initYPoints = eval(args.initYPoints)
nPts = len(initXPoints)
assert(nPts == len(initYPoints))
assert(nPts >= 2)
xyAB = numpy.zeros((nPts, 2), numpy.float64)
xyAB[:, 0] = initXPoints
xyAB[:, 1] = initYPoints
xyAB = xyAB.reshape((-1,))
print(f'number of points: {nPts}')
print(f'initial line: {xyAB}')


# integrate the line in 2D
vxyAB = numpy.zeros(xyAB.shape, numpy.float64)
def tendency(xyAB, t):
    for i in range(nPts):
        x = xyAB[i*2 + 0]
        y = xyAB[i*2 + 1]
        vxyAB[i*2 + 0] = eval(args.u_funct)
        vxyAB[i*2 + 1] = eval(args.v_funct)
    return vxyAB


# integrate the trajectories. We're creating a simple, one cell grid
# which we then advect
linePointData = vtk.vtkDoubleArray()
linePoints = vtk.vtkPoints()
lineGrid = vtk.vtkUnstructuredGrid()
lineWriter = vtk.vtkUnstructuredGridWriter()

linePointData.SetNumberOfComponents(3)
linePointData.SetNumberOfTuples(nPts)
for i in range(nPts):
    linePointData.SetTuple(i, (xyAB[i*2 + 0], xyAB[i*2 + 1], 0.0))

linePoints.SetData(linePointData)
lineGrid.SetPoints(linePoints)

abIds = vtk.vtkIdList()
abIds.SetNumberOfIds(2)
nSegs = nPts - 1
for iSeg in range(nSegs):
    i0, i1 = iSeg, iSeg + 1
    abIds.SetId(0, i0)
    abIds.SetId(1, i1)
    lineGrid.InsertNextCell(vtk.VTK_LINE, abIds)

lineWriter.SetInputData(lineGrid)
print(f'writing line to arrow_{0:05}.vtk')
lineWriter.SetFileName(f'arrow_{0:05}.vtk')
lineWriter.Update()

for iStep in range(args.numSteps):
    xyTrajectory = odeint(tendency, xyAB, [0.0, args.timeStep])
    # we're only interested in the final positions
    xyAB = xyTrajectory[1, :]
    for i in range(nPts):
        linePointData.SetTuple(i, (xyAB[i*2 + 0], xyAB[i*2 + 1], 0.0))
    print(f'writing line to arrow_{iStep:05}.vtk')
    lineWriter.SetFileName(f'arrow_{iStep:05}.vtk')
    lineWriter.Update()  





