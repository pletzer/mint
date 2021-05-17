from mint import Grid
import numpy
import sys
from pathlib import Path
import vtk


def test_create_grid():
    # create the grid
    gr = Grid()
    # 2 cells
    points = numpy.array([(0., 0., 0.),
                          (1., 0., 0.),
                          (1., 1., 0.),
                          (0., 1., 0.),
                          (1., 0., 0.),
                          (2., 0., 0.),
                          (2., 1., 0.),
                          (1., 1., 0.)]).reshape(2, 4, 3)
    gr.setPoints(points)
    gr.dump('test_create_grid.vtk')

def test_attach_data():
    # create the grid
    gr = Grid()
    # 2 cells
    points = numpy.array([(0., 0., 0.),
                          (1., 0., 0.),
                          (1., 1., 0.),
                          (0., 1., 0.),
                          (1., 0., 0.),
                          (2., 0., 0.),
                          (2., 1., 0.),
                          (1., 1., 0.)]).reshape(2, 4, 3)
    gr.setPoints(points)
    # create cell data, 3 per cell
    nDataPerCell = 3
    data = numpy.arange(0, 2*nDataPerCell, dtype=numpy.float64).reshape((2, nDataPerCell))
    gr.attach('mydata', data)
    gr.dump('test_create_grid_with_data.vtk')
    # read the data back to check the layout of the data
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName('test_create_grid_with_data.vtk')
    reader.Update()
    ugrid = reader.GetOutput()
    arr = ugrid.GetCellData().GetArray('mydata')
    assert(arr.GetNumberOfTuples() == gr.getNumberOfCells())
    assert(arr.GetNumberOfComponents() == nDataPerCell)

def test_load_grid():
    gr = Grid()
    gr.load('test_create_grid.vtk')
    ncells = gr.getNumberOfCells()
    print(f'ncells = {ncells}')
    assert(ncells == 2)

def test_load_from_ugrid_file():
    data_dir = Path(__file__).absolute().parent / '../../data'
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(data_dir / Path('cs_4.nc'))
    gr.loadFrom2DUgrid(f'{filename}:physics')
    nedges = gr.getNumberOfEdges()
    print(f'nedges = {nedges}')
    assert(nedges == 192)
    ncells = gr.getNumberOfCells()
    for icell in range(ncells):
        for iedge in range(4):
            edgeId, edgeSign = gr.getEdgeId(icell, iedge)
            nodeIds = gr.getNodeIds(icell, iedge)
            print(f'cell {icell} edge {iedge}: edgeId = {edgeId}, {edgeSign} nodeIds = {nodeIds}')
    # attaching a 3 components field to the grid
    data = numpy.array(range(ncells*4*3), numpy.float64)
    gr.attach('myData', data)

def test_edge_arc_lengths():
    data_dir = Path(__file__).absolute().parent / '../../data'
    gr = Grid()
    gr.setFlags(1, 1)
    filename = str(data_dir / Path('cs_4.nc'))
    gr.loadFrom2DUgrid(f'{filename}:physics')
    gr.computeEdgeArcLengths()
    ncells = gr.getNumberOfCells()
    for icell in range(ncells):
        for edgeIndex in range(4):
            arcLength = gr.getEdgeArcLength(icell, edgeIndex)
            print(f'cell {icell} edge {edgeIndex} edge arc length (radius=1): {arcLength}')


if __name__ == '__main__':

    test_edge_arc_lengths()
    test_attach_data()
    test_create_grid()
    test_load_grid()
    test_load_from_ugrid_file()

