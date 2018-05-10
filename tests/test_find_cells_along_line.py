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

    print 'points {} -> {} cells '.format(p0, p1),
    for i in range(cellIds.GetNumberOfIds()):
        print  '{} '.format(cellIds.GetId(i)),
    print

def test2(p0, p1, tol):

    cellIds = vtk.vtkIdList()

    ug = UgridReader('mesh_C4.nc')
    grid = ug.getUnstructuredGrid()
    loc = vtk.vtkCellLocator()
    loc.SetDataSet(grid)
    loc.BuildLocator()

    point0 = numpy.array([p0[0], p0[1], 0.0])
    point1 = numpy.array([p1[0], p1[1], 0.0])
    loc.FindCellsAlongLine(point0, point1, tol, cellIds)

    print 'points {} -> {} cells '.format(p0, p1),
    for i in range(cellIds.GetNumberOfIds()):
        print  '{} '.format(cellIds.GetId(i)),
    print

def main():
    print 'Test1' + 40*'-'
    test1([pi/2., 0.], [pi/2.+2.*pi/16., 0.], tol=0.1)
    print 'Test2' + 40*'-'
    test2([pi/2., 0.], [pi/2.+2.*pi/16., 0.], tol=1.e-10)
    print 'Test3' + 40*'-'
    test2([1.9634954084936207, 0.36548975596819283], [1.5707963267948966, 0.39269908169872414, 0.0], tol=1.e-10)

if __name__ == '__main__':
    main()

