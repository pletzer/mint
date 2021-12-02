from mint import Grid, PolylineIntegral, printLogMessages, writeLogMessages
import numpy
import pytest
import vtk


def streamFunction(p):
    # singularity is at (0, 0)
    x, y = p[:2]
    angle = numpy.arctan2(y, x)
    return angle/(2*numpy.pi)


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

            #  3--<--2
            #  |     |
            #  v     ^
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
            s0 = streamFunc(points[k, 0, :])
            s1 = streamFunc(points[k, 1, :])
            s2 = streamFunc(points[k, 2, :])
            s3 = streamFunc(points[k, 3, :])

            data[k, 0] = s1 - s0

            data[k, 1] = s2 - s1
            if data[k, 1] > 0.5:
                # points 2 and 1 are on either side of the branch
                data[k, 1] -= 1.
            elif data[k, 1] < -0.5:
                data[k, 1] += 1.

            data[k, 2] = s2 - s3

            data[k, 3] = s3 - s0
            if data[k, 3] > 0.5:
                # points 3 and 0 are on either side of the branch
                data[k, 3] -= 1.
            elif data[k, 3] < -0.5:
                data[k, 3] += 1.

            # increment the cell counter
            k += 1

    gr.setPoints(points)
    gr.dump('grid.vtk')

    return gr, data


def createLoop(xybeg, xyend, nt):

    dt = numpy.pi/float(nt)
    a = (xyend[1] - xybeg[1])/2.
    xyz = [(xybeg[0], xybeg[1], 0.), (0., xybeg[1], 0.)] + \
          [(a*numpy.cos(t), a*numpy.sin(t), 0.) \
               for t in numpy.linspace(-numpy.pi/2. + dt, numpy.pi/2. - dt, nt-1)] + \
          [(0., xyend[1], 0.), (xyend[0], xyend[1], 0.)]

    return numpy.array(xyz)


def createCircle(xycenter, nt, radius):

    nt1 = nt + 1

    xyz = numpy.zeros((nt1, 3), numpy.float64)
    ts = numpy.linspace(0., 2*numpy.pi, nt1)
    xyz[:, 0] = xycenter[0] + radius*numpy.cos(ts)
    xyz[:, 1] = xycenter[1] + radius*numpy.sin(ts)

    return xyz


def saveLineVTK(xyz, filename):

    n = xyz.shape[0]
    ptsData = vtk.vtkDoubleArray()
    ptsData.SetNumberOfComponents(3)
    ptsData.SetNumberOfTuples(n)
    ptsData.SetVoidArray(xyz, n*3, 1)

    pts = vtk.vtkPoints()
    pts.SetData(ptsData)

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions((n, 1, 1))
    sgrid.SetPoints(pts)

    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(sgrid)
    writer.Update()


def test_fluxes(nx, ny):

    # the singularity must fall inside a cell
    assert(nx % 2 == 1)
    assert(ny % 2 == 1)

    xymin = (-1., -1.)
    xymax = (1., 1.)

    dx = (xymax[0] - xymin[0])/float(nx)
    dy = (xymax[1] - xymin[1])/float(ny)
    h = min(dx, dy)

    grid, data = createGridAndData(nx=nx, ny=ny, xymin=xymin, xymax=xymax,
                                   streamFunc=streamFunction)

    results = {
        'A': {'xyz': createCircle(xycenter=(0., 0.), nt=8, radius=h),
                      'flux': float('nan'),
                      'exact': 1.0,
                     },
        'B': {'xyz': createCircle(xycenter=(0., 0.), nt=32, radius=0.7),
                       'flux': float('nan'),
                       'exact': 1.0,
                     },
        'C': {'xyz': createLoop(xybeg=(1., 1.5*h), xyend=(1., -1.5*h), nt=15),
                       'flux': float('nan'),
                       'exact': streamFunction((1., -1.5*h)) - streamFunction((1., 1.5*h)) + 1.,
                     },
        'D': {'xyz': createLoop(xybeg=(0.91, 0.55), xyend=(0.91, -0.55), nt=15),
                       'flux': float('nan'),
                       'exact': streamFunction((1., -0.55)) - streamFunction((1., 0.55)) + 1.,
                     },
    }

    for case in results:

        pli = PolylineIntegral()
        # no periodicity in x
        pli.build(grid, results[case]['xyz'], counterclock=False, periodX=0.0)

        results[case]['flux'] = pli.getIntegral(data)

        # save the contour in VTK file
        saveLineVTK(results[case]['xyz'], case + '.vtk')

        print(f'{case}: flux = {results[case]["flux"]} exact = {results[case]["exact"]} error = {results[case]["flux"] - results[case]["exact"]}')


if __name__ == '__main__':
    amax = 1.
    test_fluxes(nx=11, ny=11)
