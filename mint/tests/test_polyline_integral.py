from mint import Grid, PolylineIntegral
import numpy
import sys
from pathlib import Path


def potentialFunc(p):
    x, y = p[:2]
    return x + 2*y


def test_simple():

    # create the grid and the edge data
    gr = Grid()

    nx, ny = 3, 2
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
            data[k, 0] = potentialFunc(points[k, 1, :]) - potentialFunc(points[k, 0, :])
            data[k, 1] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 1, :])
            data[k, 2] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 3, :])
            data[k, 3] = potentialFunc(points[k, 3, :]) - potentialFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)
    gr.dump('test_polyline_integral.vtk')


    pli = PolylineIntegral()

    # create the polyline through which the flux will be integrated
    xyz = numpy.array([(0., 0., 0.),
                       (1., 0., 0.),
                       (1., 1., 0.),
                       (0., 1., 0.)])

    # no periodicity in x
    pli.build(gr, xyz, counterclock=False, periodX=0.0)

    flux = pli.getIntegral(data)
    exactFlux = potentialFunc(xyz[-1, :]) - potentialFunc(xyz[0, :])
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10

def test_partially_outside():

    # create the grid and the edge data
    gr = Grid()

    nx, ny = 3, 2
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
            data[k, 0] = potentialFunc(points[k, 1, :]) - potentialFunc(points[k, 0, :])
            data[k, 1] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 1, :])
            data[k, 2] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 3, :])
            data[k, 3] = potentialFunc(points[k, 3, :]) - potentialFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)

    pli = PolylineIntegral()

    # create the polyline through which the flux will be integrated
    xyz = numpy.array([(-0.5, 0., 0.),
                       (1., 0., 0.),
                       (1., 1., 0.),
                       (0., 1., 0.)])

    # no periodicity in x
    pli.build(gr, xyz, counterclock=False, periodX=0.0)

    flux = pli.getIntegral(data)

    # because the first point is outside the domain, only the contribution
    # stemming from the path inside the domain will be computed. Let's 
    # correct for this by moving the first point inwards
    xyz[0, 0] = 0.
    exactFlux = potentialFunc(xyz[-1, :]) - potentialFunc(xyz[0, :])
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10

def test_completely_outside():

    # create the grid and the edge data
    gr = Grid()

    nx, ny = 3, 2
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
            data[k, 0] = potentialFunc(points[k, 1, :]) - potentialFunc(points[k, 0, :])
            data[k, 1] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 1, :])
            data[k, 2] = potentialFunc(points[k, 2, :]) - potentialFunc(points[k, 3, :])
            data[k, 3] = potentialFunc(points[k, 3, :]) - potentialFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)

    pli = PolylineIntegral()

    # create the polyline through which the flux will be integrated
    xyz = numpy.array([(0., 0., 0.),
                       (-1., 0., 0.),
                       (-1., 1., 0.),
                       (0., 1., 0.)])

    # no periodicity in x
    pli.build(gr, xyz, counterclock=False, periodX=0.0)

    flux = pli.getIntegral(data)
    exactFlux = 0.0
    print(f'total flux: {flux:.3f} exact flux: {exactFlux:.3f}')
    assert abs(flux - exactFlux) < 1.e-10


if __name__ == '__main__':

    test_simple()
    test_partially_outside()
    test_completely_outside()


