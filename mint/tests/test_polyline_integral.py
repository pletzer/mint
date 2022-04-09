from errno import EBADARCH
from mint import (Grid, PolylineIntegral, printLogMessages, writeLogMessages,
                 CELL_BY_CELL_DATA, UNIQUE_EDGE_DATA)
import numpy
import pytest
from pathlib import Path


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def potentialFunc(p):
    x, y = p[:2]
    return x + 2*y


def singularPotentialFunc(p):
    x, y = p[:2]
    return numpy.arctan2(y, x)/(2.*numpy.pi)


@pytest.mark.parametrize("nx", [3])
@pytest.mark.parametrize("ny", [2])
@pytest.mark.parametrize("potFunc", [potentialFunc, singularPotentialFunc])
@pytest.mark.parametrize("xyz", [numpy.array([(1., 0., 0.),
                                              (0., 1., 0.)]),
                                 numpy.array([(0., 0., 0.),
                                              (1., 0., 0.),
                                              (1., 1., 0.),
                                              (0., 1., 0.)])])
def test_simple(nx, ny, potFunc, xyz):

    # create the grid and the edge data
    gr = Grid()

    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    data = numpy.zeros((nx*ny, 4))
    dx = 1.0 / float(nx)
    dy = 1.0 / float(ny)
    k = 0
    for i in range(nx):
        x0 = i*dx
        x1 = x0 + dx
        for j in range(ny):
            y0 = j*dy
            y1 = y0 + dy

            # node indexing
            #  3-->--2
            #  |     |
            #  ^     ^
            #  |     |
            #  0-->--1
            points[k, 0, :] = x0, y0, 0.
            points[k, 1, :] = x1, y0, 0.
            points[k, 2, :] = x1, y1, 0.
            points[k, 3, :] = x0, y1, 0.

            # edge indexing
            #     2
            #  +-->--+
            #  |     |
            # 3^     ^1
            #  |     |
            #  +-->--+
            #     0
            data[k, 0] = potFunc(points[k, 1, :]) - potFunc(points[k, 0, :])
            data[k, 1] = potFunc(points[k, 2, :]) - potFunc(points[k, 1, :])
            data[k, 2] = potFunc(points[k, 2, :]) - potFunc(points[k, 3, :])
            data[k, 3] = potFunc(points[k, 3, :]) - potFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)
    gr.dump('test_polyline_integral.vtk')

    pli = PolylineIntegral()

    pli.setGrid(gr)

    # no periodicity in x
    pli.buildLocator(numCellsPerBucket=128, periodX=0.0, enableFolding=False)

    pli.computeWeights(xyz, counterclock=False)

    flux = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)
    exactFlux = potFunc(xyz[-1, :]) - potFunc(xyz[0, :])
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10


@pytest.mark.parametrize("nx", [3])
@pytest.mark.parametrize("ny", [2])
@pytest.mark.parametrize("potFunc", [potentialFunc])
def test_partially_outside(nx, ny, potFunc):

    print('target line is partially outside the domain, expect a warning!')
    # create the grid and the edge data
    gr = Grid()

    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    data = numpy.zeros((nx*ny, 4))
    dx = 1.0 / float(nx)
    dy = 1.0 / float(ny)
    k = 0
    for i in range(nx):
        x0 = i*dx
        x1 = x0 + dx
        for j in range(ny):
            y0 = j*dy
            y1 = y0 + dy

            # node indexing
            #  3-->--2
            #  |     |
            #  ^     ^
            #  |     |
            #  0-->--1
            points[k, 0, :] = x0, y0, 0.
            points[k, 1, :] = x1, y0, 0.
            points[k, 2, :] = x1, y1, 0.
            points[k, 3, :] = x0, y1, 0.

            # edge indexing
            #     2
            #  +-->--+
            #  |     |
            # 3^     ^1
            #  |     |
            #  +-->--+
            #     0
            data[k, 0] = potFunc(points[k, 1, :]) - potFunc(points[k, 0, :])
            data[k, 1] = potFunc(points[k, 2, :]) - potFunc(points[k, 1, :])
            data[k, 2] = potFunc(points[k, 2, :]) - potFunc(points[k, 3, :])
            data[k, 3] = potFunc(points[k, 3, :]) - potFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)

    pli = PolylineIntegral()

    # create the polyline through which the flux will be integrated
    xyz = numpy.array([(-0.5, 0., 0.),
                       (1., 0., 0.),
                       (1., 1., 0.),
                       (0., 1., 0.)])

    pli.setGrid(gr)

    # no periodicity in x
    pli.buildLocator(numCellsPerBucket=128, periodX=0.0, enableFolding=False)

    pli.computeWeights(xyz, counterclock=False)

    flux = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)

    # because the first point is outside the domain, only the contribution
    # stemming from the path inside the domain will be computed. Let's
    # correct for this by moving the first point inwards
    xyz[0, 0] = 0.
    exactFlux = potentialFunc(xyz[-1, :]) - potentialFunc(xyz[0, :])
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10

    # print all the log messages 
    printLogMessages()
    # write the logs to file
    writeLogMessages("test_partially_outside_log.txt")


@pytest.mark.parametrize("nx", [3])
@pytest.mark.parametrize("ny", [2])
@pytest.mark.parametrize("potFunc", [potentialFunc])
def test_completely_outside(nx, ny, potFunc):

    print('target line is outside the domain, expect warnings!')
    # create the grid and the edge data
    gr = Grid()

    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    data = numpy.zeros((nx*ny, 4))
    dx = 1.0 / float(nx)
    dy = 1.0 / float(ny)
    k = 0
    for i in range(nx):
        x0 = i*dx
        x1 = x0 + dx
        for j in range(ny):
            y0 = j*dy
            y1 = y0 + dy

            # node indexing
            #  3-->--2
            #  |     |
            #  ^     ^
            #  |     |
            #  0-->--1
            points[k, 0, :] = x0, y0, 0.
            points[k, 1, :] = x1, y0, 0.
            points[k, 2, :] = x1, y1, 0.
            points[k, 3, :] = x0, y1, 0.

            # edge indexing
            #     2
            #  +-->--+
            #  |     |
            # 3^     ^1
            #  |     |
            #  +-->--+
            #     0
            data[k, 0] = potFunc(points[k, 1, :]) - potFunc(points[k, 0, :])
            data[k, 1] = potFunc(points[k, 2, :]) - potFunc(points[k, 1, :])
            data[k, 2] = potFunc(points[k, 2, :]) - potFunc(points[k, 3, :])
            data[k, 3] = potFunc(points[k, 3, :]) - potFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)

    pli = PolylineIntegral()
    pli.setGrid(gr)

    # create the polyline through which the flux will be integrated
    xyz = numpy.array([(0., 0., 0.),
                       (-1., 0., 0.),
                       (-1., 1., 0.),
                       (0., 1., 0.)])

    # no periodicity in x
    pli.buildLocator(numCellsPerBucket=128, periodX=0.0, enableFolding=False)

    pli.computeWeights(xyz, counterclock=False)

    flux = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)
    exactFlux = 0.0
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10


def test_identity():

    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True) # uniform lat-lon
    filename = str(DATA_DIR / Path('latlon4x2.nc'))
    meshname = 'latlon'
    grid.loadFromUgrid2DFile(f'{filename}${meshname}')
    num_edges = grid.getNumberOfEdges()
    data = numpy.array(range(0, num_edges), numpy.float64)

    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=100, periodX=0., enableFolding=False)

    xyz = numpy.array([(0., 0., 0.),
                       (90., 0., 0.),])
    pli.computeWeights(xyz, counterclock=False)
    flux = pli.getIntegral(data=data, placement=UNIQUE_EDGE_DATA)

    print(f'flux = {flux}')
    # assert abs(flux - .) < 1.e-10


if __name__ == '__main__':

    test_identity()

    # # polyline through which the line integral will be computed
    # xyz = numpy.array([(1., 0., 0.),
    #                    (0., 1., 0.)])
    # test_simple(3, 2, singularPotentialFunc, xyz)

    # xyz = numpy.array([(0., 0., 0.),
    #                    (1., 0., 0.),
    #                    (1., 1., 0.),
    #                    (0., 1., 0.)])
    # test_simple(3, 2, potentialFunc, xyz)

    # test_partially_outside(2, 3, potentialFunc)
    # test_completely_outside(2, 3, potentialFunc)
