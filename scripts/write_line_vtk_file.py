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

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(args.o)
writer.SetInputData(grid)
writer.Update()
