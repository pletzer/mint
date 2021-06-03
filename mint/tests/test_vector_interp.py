from mint import Grid
from mint import VectorInterp
import numpy
import sys
from pathlib import Path
import vtk


def test_one_cell():
    # create the grid, a single slanted cell
    gr = Grid()
    points = numpy.array([(0., 0., 0.),
                          (1., 0., 0.),
                          (1., 1., 0.),
                          (0.8, 1., 0.)]).reshape(1, 4, 3)
    gr.setPoints(points)

    # create the interpolator
    vi = VectorInterp()
    vi.setGrid(gr)
    vi.buildLocator(numCellsPerBucket=10, periodX=0.)

    # target points
    x = numpy.linspace(0., 1., 11)
    y = numpy.linspace(0., 1., 21)
    xx, yy = numpy.meshgrid(x, y)
    xflat, yflat = xx.flat, yy.flat
    targetPoints = numpy.array([(xflat[i], yflat[i], 0.0) for i in range(len(xflat))])
    numBad = vi.findPoints(targetPoints, tol2=1.e-10)
    assert(numBad == 92)

    for data in (1., 0., 0., 0.), (0., 1., 0., 0.), (0., 0., 1., 0.), (0., 0., 0., 1.):
        vectors = vi.getVectors(numpy.array(data))
        # NEED TO ADD SOME CHECKS HERE


if __name__ == '__main__':

    test_one_cell()

