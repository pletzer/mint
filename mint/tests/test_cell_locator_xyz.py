"""
Regression test for a bug in vmtLonLatCellLocator::getBucketId (formerly
vmtCellLocator, see vmtCellLocator.h/vmtLonLatCellLocator.h): it bins a
point into the locator's spatial hash using ONLY its first two (x, y)
coordinates --

    x[0] = (point[0] - xmin[0]) / (xmax[0] - xmin[0])
    x[1] = (point[1] - xmin[1]) / (xmax[1] - xmin[1])
    ... m, n derived from x[0], x[1] only, point[2] (z) is never read ...

-- the same "2D coordinate-chart" assumption pattern already found and
fixed in mnt_vectorinterp__getTangentVectors (mntVectorInterp.h) and
isPointInQuad (vmtLonLatCellLocator.cpp), but in a THIRD, separate place,
and a more fundamental one: the bucket assignment is a PRE-FILTER that
runs before any per-cell containment test, so if a point's true
containing cell isn't registered in the bucket the point's own (x, y)
falls into, FindCell never even considers that cell -- no amount of
fixing the downstream (already-3D-aware) containment/interpolation math
can recover from that.

FIX: rather than trying to make vmtLonLatCellLocator's own bucket grid
3D-aware, this class's periodicity/pole-folding machinery is genuinely
unneeded for a plain embedded (x, y, z) mesh -- so VectorInterp.buildLocator
now takes a useXYZLocator flag (default False, preserving the exact prior
behaviour) that, when True, builds a NEW class, vmtXYZCellLocator, instead:
a thin wrapper around a standard vtkStaticCellLocator, which indexes on all
3 coordinates and has no periodicity assumptions to get wrong. This test
passes by opting into that (see buildLocator's useXYZLocator=True call
below).

Mechanism: BuildLocator registers each cell into the bucket(s) of its OWN
4 corners (not its full interior). For a genuinely 3D-embedded, CURVED
surface, a cell can be "edge-on" to the (x, y) plane -- large z-extent,
small/awkward (x, y) footprint -- so a target point that is a bona fide
interior (bilinear, convex-combination) point of the cell can round to a
DIFFERENT (x, y) bucket than any of that cell's own corners occupy,
especially once the panel is subdivided finely enough (M=64) that
floating-point rounding near a bucket boundary tips the balance. This is
exactly what the module docstring's own "WARNING: this could fail if the
buckets are much smaller than some cells!" refers to -- except here it is
the projection, not the cell/bucket size ratio, that is the real cause.

Found via scripts/xyz_pole_panel_validation.py's convergence study on a
genuine (x, y, z) unit-sphere '+X' panel at M=64: 2 of 4096 target points
failed findPoints even though they are ordinary interior points of
ordinary (non-degenerate) cells -- confirmed in isolation (a single-cell
grid built from just the failing cell's own 4 corners locates the same
target point without any problem; only the FULL 4096-cell grid's bucket
structure reproduces the failure -- smaller sub-grids around the same
cell, with a coarser bucket structure, do not).

This test was written to fail first (mirroring the same
write-failing-test-first pattern used for the mntVectorInterp.h /
vmtLonLatCellLocator.cpp fixes in test_vector_interp_xyz.py), confirmed to
fail against vmtLonLatCellLocator, and now passes via vmtXYZCellLocator.
"""
import numpy

from mint import Grid, VectorInterp


def buildPlusXPanelXYZ(M):
    """
    M x M cells covering the '+X' face of a cube, gnomonically projected
    onto the unit sphere -- genuine embedded (x, y, z) points, no lon/lat
    at all (same construction as scripts/xyz_pole_panel_validation.py's
    buildPlusZPanelXYZ, just centred on '+X' instead of '+Z'; duplicated
    here to keep this test file self-contained).
    """
    normal, u_hat, v_hat = numpy.array([1., 0., 0.]), numpy.array([0., 1., 0.]), numpy.array([0., 0., 1.])
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = numpy.meshgrid(edges, edges, indexing='ij')
    cube_pts = (normal[None, None, :] + U[..., None] * u_hat[None, None, :]
                + V[..., None] * v_hat[None, None, :])
    norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
    sphere_pts = cube_pts / norm

    points = numpy.empty((M * M, 4, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            points[icell, 0] = sphere_pts[i, j]
            points[icell, 1] = sphere_pts[i + 1, j]
            points[icell, 2] = sphere_pts[i + 1, j + 1]
            points[icell, 3] = sphere_pts[i, j + 1]
            icell += 1
    return points


def buildXYZGrid(points):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=0)
    grid.setPoints(points)
    return grid


def test_locator_finds_interior_point_on_curved_xyz_panel_M64():
    """
    Two known-bad (icell, target) pairs on the M=64 '+X' xyz panel -- both
    targets are the bilinear (0.25, 0.75) point of an ordinary interior
    cell (not a boundary/pole cell), yet findPoints reports them as
    outside the grid.
    """
    M = 64
    points = buildPlusXPanelXYZ(M)

    grid = buildXYZGrid(points)
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=128, periodX=0.0, enableFolding=False, useXYZLocator=True)

    # (icell, bilinear xi/eta) -- targets computed and renormalized onto
    # the unit sphere the same way xyz_pole_panel_validation.py's
    # bilinearCellTargets does
    xi, eta = 0.25, 0.75
    knownBadCells = [2580, 2636]

    for icell in knownBadCells:
        p0, p1, p2, p3 = points[icell]
        target = ((1 - xi) * (1 - eta) * p0 + xi * (1 - eta) * p1
                  + xi * eta * p2 + (1 - xi) * eta * p3)
        target = target / numpy.linalg.norm(target)

        numBad = vi.findPoints(numpy.array([target]), tol2=1.e-6)
        assert numBad == 0, (
            f'cell {icell}: interior target point {target} was not found by the locator '
            f'-- see this module\'s docstring (vmtCellLocator::getBucketId only bins on x,y)')


if __name__ == '__main__':
    test_locator_finds_interior_point_on_curved_xyz_panel_M64()
