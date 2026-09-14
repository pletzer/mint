"""
Regression tests for mnt_vectorinterp__getTangentVectors (mntVectorInterp.h)
and isPointInQuad (vmtCellLocator.cpp) no longer assuming z=0 -- i.e. that
mint.VectorInterp/mint.Grid work correctly on a genuinely 3D-embedded 2D
mesh (a cell whose 4 corners do NOT lie in a z=const plane), not just on
the (lon, lat, elev=0) grids every other test in this repo uses.

Before the fix, both functions computed their Jacobian/normal from only the
first two (x, y) coordinates of each point (crossDotZHat / dx*dy-dy*dx),
implicitly assuming normal=(0,0,1). For a cell lying in a plane where that's
badly wrong -- e.g. entirely in the x=0 (y-z) plane, where every edge has
zero x-component -- the OLD code would compute a Jacobian of EXACTLY ZERO
for every corner (crossDotZHat(a, b) = a[0]*b[1] - a[1]*b[0] is identically
0 when a[0]=b[0]=0), triggering mnt_vectorinterp__getTangentVectors' "bad
cell" warning and producing NaN/garbage (division by zero) in the
reconstructed vectors -- for a perfectly valid, non-degenerate quad. That
is the sharpest, most direct illustration of the bug this file guards
against, so test_x_const_quad_is_not_degenerate targets it directly.

See also test_vector_interp.py's test_accuracy(), which independently
reimplements (in pure Python, not calling into mint at all) the same
normal/Jacobian formula this fix ports into the real C++ library, and
already demonstrated (via that separate Python model) that Cartesian
coordinates should outperform lon-lat -- this file is what makes the
actual mint.VectorInterp C++ code live up to that, rather than a parallel
Python implementation proving it in principle.
"""
import numpy

from mint import Grid, VectorInterp, CELL_BY_CELL_DATA


def buildSingleCellGrid(p0, p1, p2, p3):
    """
    One-cell mint.Grid from 4 corners (each a length-3 array/tuple), in the

        3-->--2
        |     |
        ^     ^
        |     |
        0-->--1

    convention used throughout this codebase.
    """
    points = numpy.array([[p0, p1, p2, p3]], dtype=numpy.float64)
    grid = Grid()
    grid.setPoints(points)
    return grid


def buildInterpolator(grid):
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=1, periodX=0., enableFolding=False, useXYZLocator=True)
    return vi


def constantFieldEdgeData(p0, p1, p2, p3, v):
    """
    Cell-by-cell edge circulations for a CONSTANT vector field v, in the
    "conceptual xi/eta-increasing direction" convention
    mnt_vectorinterp__getEdgeVectorsFromCellByCellData actually expects:
    edges 0 and 2 both increasing xi (corner0->corner1 and corner3->corner2
    senses), edges 1 and 3 both increasing eta (corner1->corner2 and
    corner0->corner3 senses) -- NOT the literal corner_i->corner_{i+1}
    traversal for i0=2,3, which runs backwards relative to that. Same
    `sign = 1 - 2*(i0 // 2)` convention as
    vector_potential_reconstruction_convergence.buildExactEdgeData.

    A constant field is exactly reproducible by the lowest-order (k=1)
    Whitney/Nedelec element on any affine (parallelogram) cell, at every
    interior point -- not just the centroid -- so this is a strong,
    tolerance-tight check with a trivially known right answer, regardless
    of how the cell is oriented in 3D.
    """
    v = numpy.asarray(v, dtype=numpy.float64)
    corners = [numpy.asarray(p, dtype=numpy.float64) for p in (p0, p1, p2, p3)]
    data = numpy.zeros((1, 4), dtype=numpy.float64)
    for i0 in range(4):
        i1 = (i0 + 1) % 4
        sign = 1 - 2 * (i0 // 2)
        data[0, i0] = sign * numpy.dot(v, corners[i1] - corners[i0])
    return data


TEST_XI_ETAS = [(0.5, 0.5), (0.2, 0.7), (0.9, 0.1), (0.5, 0.1), (0.1, 0.5)]


def assertConstantFieldReproducedExactly(grid, data, v, tol=1.e-10):
    vi = buildInterpolator(grid)
    for xi, eta in TEST_XI_ETAS:
        # bilinear position, matching the grid's own corner ordering
        pts = grid.getPoints()[0]
        target = ((1 - xi) * (1 - eta) * pts[0] + xi * (1 - eta) * pts[1]
                  + xi * eta * pts[2] + (1 - xi) * eta * pts[3])
        numBad = vi.findPoints(numpy.array([target]), tol2=1.e-10)
        assert numBad == 0, f'target at (xi, eta)=({xi}, {eta}) not found in its own cell'

        vec = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)[0]
        assert numpy.all(numpy.isfinite(vec)), (
            f'non-finite reconstruction at (xi, eta)=({xi}, {eta}): {vec} '
            f'-- likely a zero/degenerate Jacobian (the bug this test guards against)')
        err = numpy.linalg.norm(vec - numpy.asarray(v, dtype=numpy.float64))
        assert err < tol, f'(xi, eta)=({xi}, {eta}): expected {v}, got {vec}, err={err:.3e}'


def test_x_const_quad_is_not_degenerate():
    """
    A unit square lying entirely in the x=0 (y-z) plane -- every edge has
    zero x-component, so the OLD (x,y)-only Jacobian (crossDotZHat) was
    IDENTICALLY ZERO here for every corner, regardless of cell size or
    shape. This is the sharpest possible illustration of the bug: not a
    convergence-order degradation, but an outright zero/NaN.
    """
    p0, p1, p2, p3 = (0., 0., 0.), (0., 1., 0.), (0., 1., 1.), (0., 0., 1.)
    v = (0., 1., 0.)  # uniform, in the plane of the quad
    grid = buildSingleCellGrid(p0, p1, p2, p3)
    data = constantFieldEdgeData(p0, p1, p2, p3, v)
    assertConstantFieldReproducedExactly(grid, data, v)


def test_y_const_quad_is_not_degenerate():
    """Same idea, a unit square in the y=0 (x-z) plane."""
    p0, p1, p2, p3 = (0., 0., 0.), (1., 0., 0.), (1., 0., 1.), (0., 0., 1.)
    v = (1., 0., 0.)
    grid = buildSingleCellGrid(p0, p1, p2, p3)
    data = constantFieldEdgeData(p0, p1, p2, p3, v)
    assertConstantFieldReproducedExactly(grid, data, v)


def test_tilted_quad_arbitrary_normal():
    """
    A unit square in a plane whose normal, (1,1,1)/sqrt(3), is not aligned
    with any single coordinate axis -- unlike the x=0/y=0 cases above,
    this doesn't just swap which axis gets "dropped", it checks the
    general (Newell's-method / cross-product) normal computation itself.
    """
    n = numpy.array([1., 1., 1.]) / numpy.sqrt(3.)
    e1 = numpy.array([1., -1., 0.]) / numpy.sqrt(2.)
    e2 = numpy.cross(n, e1)
    assert abs(numpy.linalg.norm(e2) - 1.) < 1.e-12  # sanity check on the construction itself

    base = numpy.array([1., 1., 1.])
    p0, p1, p2, p3 = base, base + e1, base + e1 + e2, base + e2
    v = e1 + 0.5 * e2  # an arbitrary constant vector, in the plane of the quad

    grid = buildSingleCellGrid(p0, p1, p2, p3)
    data = constantFieldEdgeData(p0, p1, p2, p3, v)
    assertConstantFieldReproducedExactly(grid, data, v)


def test_isPointInQuad_rejects_point_outside_x_const_quad():
    """
    Companion check for the isPointInQuad fix specifically (not just the
    VectorInterp Jacobian): a point far outside the x=0 quad's own (y, z)
    footprint should NOT be found in it, even though it's numerically
    "flat" against the quad in x (x=0 for both the quad and, deliberately,
    the test point) -- the old (x,y)-only test could get this right or
    wrong somewhat by accident depending on which 2 axes it kept; this
    checks the fixed, normal-aware version directly.
    """
    p0, p1, p2, p3 = (0., 0., 0.), (0., 1., 0.), (0., 1., 1.), (0., 0., 1.)
    grid = buildSingleCellGrid(p0, p1, p2, p3)
    vi = buildInterpolator(grid)

    outside = numpy.array([[0., 5., 5.]])  # same x=0 plane, well outside the unit square
    numBad = vi.findPoints(outside, tol2=1.e-10)
    assert numBad == 1, 'point well outside the quad (in y, z) was incorrectly found inside it'

    inside = numpy.array([[0., 0.5, 0.5]])
    numBad = vi.findPoints(inside, tol2=1.e-10)
    assert numBad == 0, 'point at the quad centre was not found inside it'


if __name__ == '__main__':
    test_x_const_quad_is_not_degenerate()
    test_y_const_quad_is_not_degenerate()
    test_tilted_quad_arbitrary_normal()
    test_isPointInQuad_rejects_point_outside_x_const_quad()
