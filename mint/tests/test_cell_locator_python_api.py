"""
Regression test for the Python-level cell locator API (mint.LonLatCellLocator,
mint.XYZCellLocator, VectorInterp.setLocator) added on top of the C API in
src/mntVmtCellLocator.h -- before this, a Python caller could only get a
locator via VectorInterp.buildLocator, which builds AND owns one itself;
these let you build and configure one yourself (e.g. to tune it beyond what
buildLocator exposes, or to share one locator across several VectorInterp
objects) and hand it over with VectorInterp.setLocator, which takes it as a
borrowed reference (VectorInterp does not delete it).
"""
import numpy

from mint import Grid, VectorInterp, LonLatCellLocator, XYZCellLocator, CELL_BY_CELL_DATA


def test_xyz_locator_via_setLocator():
    """
    An XYZCellLocator built and configured directly (not via
    VectorInterp.buildLocator), on a cell lying entirely in the x=0 plane --
    see test_vector_interp_xyz.py's test_x_const_quad_is_not_degenerate for
    why that's a meaningful, not arbitrary, choice of cell.
    """
    p0, p1, p2, p3 = (0., 0., 0.), (0., 1., 0.), (0., 1., 1.), (0., 0., 1.)
    points = numpy.array([[p0, p1, p2, p3]])

    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=0)
    grid.setPoints(points)

    loc = XYZCellLocator()
    loc.setDataSet(grid)
    loc.setNumberOfCellsPerBucket(1)
    loc.buildLocator()

    vi = VectorInterp()
    vi.setGrid(grid)
    vi.setLocator(loc)

    target = numpy.array([[0., 0.5, 0.5]])
    numBad = vi.findPoints(target, tol2=1.e-10)
    assert numBad == 0

    v = numpy.array([0., 1., 0.])
    corners = [numpy.array(p) for p in (p0, p1, p2, p3)]
    data = numpy.zeros((1, 4))
    for i0 in range(4):
        i1 = (i0 + 1) % 4
        sign = 1 - 2 * (i0 // 2)
        data[0, i0] = sign * numpy.dot(v, corners[i1] - corners[i0])

    vec = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)[0]
    assert numpy.linalg.norm(vec - v) < 1.e-10


def test_lonlat_locator_via_setLocator():
    """A LonLatCellLocator built and configured directly, on a plain flat quad."""
    v0, v1, v2, v3 = (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.)
    points = numpy.array([[v0, v1, v2, v3]])

    grid = Grid()
    grid.setPoints(points)

    loc = LonLatCellLocator()
    loc.setDataSet(grid)
    loc.setNumberOfCellsPerBucket(1)
    loc.setPeriodicityLengthX(0.)
    loc.buildLocator()

    vi = VectorInterp()
    vi.setGrid(grid)
    vi.setLocator(loc)

    numBad = vi.findPoints(numpy.array([[0.5, 0.5, 0.]]), tol2=1.e-10)
    assert numBad == 0


if __name__ == '__main__':
    test_xyz_locator_via_setLocator()
    test_lonlat_locator_via_setLocator()
