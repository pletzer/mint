from ugrid_reader import UgridReader
import vtk
from math import pi
import numpy

def test1(p0, p1, tol):

    cellIds = vtk.vtkIdList()

    ug = UgridReader('mesh_C4.nc')
    grid = ug.getUnstructuredGrid()
    loc = vtk.vtkCellLocator()
    loc.SetDataSet(grid)
    loc.BuildLocator()

    point0 = numpy.array([p0[0], p0[1], 0.0])
    point1 = numpy.array([p1[0], p1[1], 0.0])
    loc.FindCellsAlongLine(point0, point1, tol, cellIds)

    for i in range(cellIds.GetNumberOfIds()):
        print 'points {} -> {} cell {}'.format(p0, p1, cellIds.GetId(i))

def main():
    test1([pi/2., 0.], [pi/2.+2.*pi/16., 0.], tol=1.e-2)

if __name__ == '__main__':
    main()

