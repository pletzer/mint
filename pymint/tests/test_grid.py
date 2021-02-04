from mintregrid import Grid
import numpy
import sys
from pathlib import Path


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


def test_load_grid():

    gr = Grid()

    gr.load('test_create_grid.vtk')

    ncells = gr.getNumberOfCells()
    print(f'ncells = {ncells}')
    assert ncells == 2


def test_load_from_ugrid_file(data_dir):

    gr = Grid()

    gr.setFlags(1, 1)

    filename = str(data_dir / Path('cs_4.nc'))
    gr.loadFrom2DUgrid(f'{filename}:physics')

    nedges = gr.getNumberOfEdges()
    print(f'nedges = {nedges}')
    assert nedges == 192

    ncells = gr.getNumberOfCells()
    for icell in range(ncells):
        for iedge in range(4):
            edgeId, edgeSign = gr.getEdgeId(icell, iedge)
            nodeIds = gr.getNodeIds(icell, iedge)
            print(f'cell {icell} edge {iedge}: edgeId = {edgeId}, {edgeSign} nodeIds = {nodeIds}')

    # attaching a 3 components field to the grid
    data = numpy.array(range(ncells*4*3), numpy.float64)
    gr.attach('myData', data)


if __name__ == '__main__':

    data_dir = Path('./data')
    if len(sys.argv) >= 2:
        data_dir = Path(sys.argv[1])

    test_create_grid()
    test_load_grid()
    test_load_from_ugrid_file(data_dir)


