import vtk
import argparse
import numpy

parser = argparse.ArgumentParser(description='Write line point data to file')
parser.add_argument('-p', type=str, default="(0., 0.),(1., 0.)", 
                    help='Interlaced xy points')
parser.add_argument('-o', type=str, default="line.vtk", 
                    help='Output file')
args = parser.parse_args()

xy = numpy.array(eval('[' + args.p + ']'))
npts = xy.shape[0]
ncells = npts - 1
xyz = numpy.zeros((npts, 3), numpy.float64)
xyz[:, 0] = xy[:, 0]
xyz[:, 1] = xy[:, 1]

pointData = vtk.vtkDoubleArray()
pointData.SetNumberOfComponents(3)
pointData.SetNumberOfTuples(npts)
pointData.SetVoidArray(xyz, 3*npts, 1)

points = vtk.vtkPoints()
points.SetNumberOfPoints(npts)
points.SetData(pointData)

grid = vtk.vtkUnstructuredGrid()
grid.SetPoints(points)
grid.Allocate(ncells, 1)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(2)
for i in range(ncells):
	ptIds.SetId(0, i)
	ptIds.SetId(1, i + 1)
	grid.InsertNextCell(vtk.VTK_LINE, ptIds)

cellLocator = vtk.vtkCellLocator()
cellLocator.SetDataSet(grid)
cellLocator.BuildLocator()
cell = vtk.vtkGenericCell()
subId = vtk.mutable(0)
pcoords = numpy.zeros((3,), numpy.float64)
weights = numpy.zeros((4,), numpy.float64)
tol2 = 1.e-12
twopi = 2*numpy.pi
ptIds = vtk.vtkIdList()

def vectorField(xy):
	# find the cell and the pcoords
	xyz = numpy.array([xy[0], xy[1], 0.])
	cellId = cellLocator.FindCell(xyz, tol2, cell, pcoords, weights)
	if cellId < 0: print 'ERROR cellId = ', cellId

	grid.GetCellPoints(cellId, ptIds)

	# interpolate the vector field at that position
	vec = numpy.zeros((3,), numpy.float64)
	vertex = numpy.zeros((3,), numpy.float64)
	# iterate over the cell vertices
	for i in range(ptIds.GetNumberOfIds()):
		ptId = ptIds.GetId(i)
		grid.GetPoint(ptId, vertex)
		x, y, z = vertex
		r2 = x**2 + y**2
		v = numpy.array([x, y, 0.])/(twopi*r2)
		vec += weights[i] * v
	return vec

def integratedFlux(xy0, xy1):
	xyMid = 0.5*(xy0 + xy1)
	vec = vectorField(xyMid)
	ds = numpy.array([xy1[1] - xy0[1], -(xy1[0] - xy0[0]), 0.0])
	return numpy.dot(vec, ds)

flux = 0.0
for i in range(ncells):
	xy0 = xy[i + 0, :]
	xy1 = xy[i + 1, :]
	flux += integratedFlux(xy0, xy1)

print 'flux = {:16.12f}'.format(flux)
