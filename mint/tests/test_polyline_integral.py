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

    xyz = numpy.array([(0., 90., 0.),
                       (360., 90., 0.),])
    pli.computeWeights(xyz, counterclock=False)
    flux = pli.getIntegral(data=data, placement=UNIQUE_EDGE_DATA)

    print(f'flux = {flux}')
    assert abs(flux - 38.) < 1.e-10


def _lonlat2xyz(lon, lat):
    lon, lat = numpy.radians(lon), numpy.radians(lat)
    return numpy.stack([numpy.cos(lat) * numpy.cos(lon),
                        numpy.cos(lat) * numpy.sin(lon),
                        numpy.sin(lat)], -1)


def test_segment_across_periodic_boundary():
    """
    A segment crossing the grid's minimum longitude must give the same
    integral as the same segment shifted by 360 deg (issue #129). After
    fixLonAcrossDateline the FV3 C24 grid spans lon [-10, 440].
    """
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1, degrees=True)
    grid.loadFromFV32DFile(str(DATA_DIR / Path('C24_SCRIP_desc.181018.nc')))
    pts = grid.getPoints()
    assert abs(pts[..., 0].min() + 10.) < 1.e-8

    # midpoint-rule edge integrals of v = omega x r, a field with non-zero
    # circulation (a gradient field would hide the error)
    omega = numpy.array([0.3, 0.5, 0.8])
    x = _lonlat2xyz(pts[..., 0], pts[..., 1])
    a, b = x[:, [0, 1, 3, 0]], x[:, [1, 2, 2, 3]]
    mid = (a + b) / numpy.linalg.norm(a + b, axis=-1, keepdims=True)
    data = numpy.sum(numpy.cross(omega, mid) * (b - a), -1).ravel()

    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=100, periodX=360., enableFolding=True)

    def integral(points):
        xyz = numpy.zeros((len(points), 3))
        xyz[:, :2] = points
        pli.computeWeights(xyz, counterclock=False)
        return pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)

    # open segment crossing lon = -10
    flux = integral([(-12., 30.), (0., 30.)])
    flux_shifted = integral([(348., 30.), (360., 30.)])
    print(f'flux = {flux} shifted = {flux_shifted}')
    assert abs(flux - flux_shifted) < 1.e-10
    # exact value 0.1034; the bug added a full latitude circle, -3.77
    assert abs(flux - 0.1034) < 1.e-3

    # closed path crossing lon = -10
    gam = numpy.linspace(0., 2*numpy.pi, 21)
    r = numpy.degrees(numpy.arcsin(1./numpy.sqrt(3.)))
    circle = numpy.stack([r*numpy.sin(gam), r*numpy.cos(gam)], -1)
    flux = integral(circle + (-10., 0.))
    flux_shifted = integral(circle + (350., 0.))
    print(f'closed path flux = {flux} shifted = {flux_shifted}')
    assert abs(flux - flux_shifted) < 1.e-10


if __name__ == '__main__':

    test_identity()

    # polyline through which the line integral will be computed
    xyz = numpy.array([(1., 0., 0.),
                       (0., 1., 0.)])
    test_simple(3, 2, singularPotentialFunc, xyz)

    xyz = numpy.array([(0., 0., 0.),
                       (1., 0., 0.),
                       (1., 1., 0.),
                       (0., 1., 0.)])
    test_simple(3, 2, potentialFunc, xyz)

    test_partially_outside(2, 3, potentialFunc)
    test_completely_outside(2, 3, potentialFunc)
