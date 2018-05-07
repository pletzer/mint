from ugrid_reader import UgridReader
import vtk
from math import pi
import numpy

def test1(p0, p1, tol):

    eps = 1.73654365e-12

    cellIds = vtk.vtkIdList()

    ug = UgridReader('mesh_C4.nc')
    grid = ug.getUnstructuredGrid()
    loc = vtk.vtkCellLocator()
    loc.SetDataSet(grid)
    loc.BuildLocator()

    point0 = numpy.array([p0[0] - eps * 1.86512432134, p0[1] + eps * 2.76354653243, 0.0])
    point1 = numpy.array([p1[0] + eps * 1.96524543545, p1[1] - eps * 0.82875646565, 0.0])
    loc.FindCellsAlongLine(point0, point1, tol, cellIds)

    for i in range(cellIds.GetNumberOfIds()):
        print 'points {} -> {} cell {}'.format(p0, p1, cellIds.GetId(i))

def test2(p0, p1, tol):

    eps = 0.0

    cellIds = vtk.vtkIdList()

    ug = UgridReader('mesh_C4.nc')
    grid = ug.getUnstructuredGrid()
    loc = vtk.vtkCellLocator()
    loc.SetDataSet(grid)
    loc.BuildLocator()

    point0 = numpy.array([p0[0] - eps * 1.86512432134, p0[1] + eps * 2.76354653243, 0.0])
    point1 = numpy.array([p1[0] + eps * 1.96524543545, p1[1] - eps * 0.82875646565, 0.0])
    loc.FindCellsAlongLine(point0, point1, tol, cellIds)

    for i in range(cellIds.GetNumberOfIds()):
        print 'points {} -> {} cell {}'.format(p0, p1, cellIds.GetId(i))

def main():
    print 'Test1' + 40*'-'
    test1([pi/2., 0.], [pi/2.+2.*pi/16., 0.], tol=0.1)
    print 'Test2' + 40*'-'
    test2([pi/2., 0.], [pi/2.+2.*pi/16., 0.], tol=0.1)

if __name__ == '__main__':
    main()

