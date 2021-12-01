from mint import Grid, PolylineIntegral, printLogMessages, writeLogMessages
import numpy
import pytest
from matplotlib import pylab


def streamFunction(p):
    x, y = p[:2]
    return numpy.arctan2(y, x)/(numpy.pi/2.)


def createGridAndData(nx, ny, xymin, xymax, streamFunc):

    # make sure the singularity falls inside a cell since the 
    # streamfunction would be ill defined there
    # assert(nx % 2 != 0)
    # assert(ny % 2 != 0)

    gr = Grid()

    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    data = numpy.zeros((nx*ny, 4))
    dx = (xymax[0] - xymin[0]) / float(nx)
    dy = (xymax[1] - xymin[1]) / float(ny)
    k = 0
    for i in range(nx):
        x0 = xymin[0] + i*dx
        x1 = x0 + dx
        for j in range(ny):
            y0 = xymin[1] + j*dy
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
            data[k, 0] = streamFunc(points[k, 1, :]) - streamFunc(points[k, 0, :])
            data[k, 1] = streamFunc(points[k, 2, :]) - streamFunc(points[k, 1, :])
            data[k, 2] = streamFunc(points[k, 2, :]) - streamFunc(points[k, 3, :])
            data[k, 3] = streamFunc(points[k, 3, :]) - streamFunc(points[k, 0, :])

            # increment the cell counter
            k += 1

    gr.setPoints(points)

    return gr, data


def createTarget(xycenter, nt, radius):

    nt1 = nt + 1

    xyz = numpy.zeros((nt1, 3), numpy.float64)
    ts = numpy.linspace(0., numpy.pi/2., nt1)
    xyz[:, 0] = xycenter[0] + radius*numpy.cos(ts)
    xyz[:, 1] = xycenter[1] + radius*numpy.sin(ts)

    return xyz


def test_fluxes(nx, ny, xymin, xymax, radius, nt, ncontours, plot=False):

    grid, data = createGridAndData(nx, ny, xymin, xymax, streamFunction)
    dx = (xymax[0] - xymin[0])/float(nx)
    dy = (xymax[1] - xymin[1])/float(ny)
    h = min(dx, dy)

    # if plot:
    #     for i in range(nx + 1):
    #         pylab.plot([i*dx, i*dx], [xymin[1], xymax[1]], 'k--')
    #     for j in range(ny + 1):
    #         pylab.plot([xymin[0], xymax[0]], [j*dy, j*dy], 'k--')
    # pylab.show()

    radii = radius * numpy.linspace(1./float(ncontours), 1., ncontours)
    fluxes = numpy.zeros((ncontours,), numpy.float64)
    for i in range(ncontours):

        xycenter = (0.0, 0.0)
        xyz = createTarget(xycenter, nt, radii[i])

        # if plot:
        #     pylab.plot(xyz[:, 0], xyz[:, 1], 'g-')
    
        pli = PolylineIntegral()
        # no periodicity in x
        pli.build(grid, xyz, counterclock=False, periodX=0.0)

        fluxes[i] = pli.getIntegral(data)

    if plot:

        # pylab.axes().set_aspect('equal', 'datalim')
        # pylab.show()

        # plot the flux and its error
        pylab.plot(radii/h, fluxes, 'b-')
        error = fluxes - 1.
        print(f'error = {error}')
        pylab.plot(radii/h, error, 'r--')
        pylab.xlabel('a/h')
        pylab.legend(['computed flux', 'numerical error'])
        pylab.plot(radii/h, fluxes, 'ko', radii/h, error, 'kx', markerfacecolor='none')
        pylab.show()


if __name__ == '__main__':
    eps = 1.23e-15
    amax = 4.
    test_fluxes(nx=4, ny=4, xymin=(0.+eps, 0.+eps), xymax=(amax, amax), radius=amax, nt=32, ncontours=40, plot=True)
