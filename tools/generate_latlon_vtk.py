import vtk
import numpy
import sys
import argparse

from math import sin, cos, pi

"""
Generate grid and edge data on uniform grid and save result in UGRID file
"""

parser = argparse.ArgumentParser(description='Generate the output file storing the lat-lon grid and edge data in VTK unstructured grid format')
parser.add_argument('-o', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-nx', default=1, type=int, 
                    help='Number of longitude cells')
parser.add_argument('-ny', default=1, type=int, 
                    help='Number of latitude cells')
parser.add_argument('-s', type=str, dest='stream_funct', default='sin(pi*x/180.)*cos(pi*y/180.)', help='Stream function of x (longitude in deg) and y (latitude in deg) used for setting the edge integrals')
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
edata = vtk.vtkDoubleArray()
writer = vtk.vtkUnstructuredGridWriter()

pdata.SetNumberOfComponents(3)
pdata.SetNumberOfTuples(npoints)

sdata.SetNumberOfComponents(1)
sdata.SetNumberOfTuples(npoints)
sdata.SetName('stream_function')

edata.SetNumberOfComponents(4)
edata.SetNumberOfTuples(nfaces)
edata.SetName('edge_velocity')


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
		edata.SetTuple(iface, (s1 - s0, s2 - s1, s2 - s3, s3 - s0))
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
grid.GetCellData().AddArray(edata)

# write to file
writer.SetFileName(args.grid_file)
writer.SetInputData(grid)
writer.Update()


