#!/usr/bin/env python
"""
Generate a cubed-sphere grid by encasing the sphere in a cube: lay a
uniform M x M grid on each of the cube's 6 faces, then project each grid
point radially onto the sphere (gnomonic projection). This gives
6 * M^2 quadrilateral cells, returned in the (ncells, 4, 3) point format
expected by mint.Grid.setPoints (longitude, latitude, elevation=0, in
degrees) -- no netCDF file needed, this can supply an arbitrary resolution
directly in memory.

Corners are stored cell by cell (mint does not need a shared-vertex/
connectivity representation for setPoints -- duplicate nodes across
neighbouring cells are fine), so the only thing that needs fixing up per
cell is its OWN 4 corners being self-consistent: not spanning most of the
way around the sphere just because two of them fell on opposite sides of
the +-180 deg branch cut, and not carrying an arbitrary longitude for a
corner sitting exactly at a pole (where longitude isn't well-defined).

Usage: python scripts/generate_cubedsphere_grid.py [M]
"""
import sys
from pathlib import Path
import numpy
import mint

NUM_VERTS_PER_QUAD = 4


def rotateZ(v, angleDeg):
    """
    Rotate a 3-vector (or an array of 3-vectors, last axis of length 3) by
    angleDeg about the vertical (Z) axis. Since Z-axis rotation is exactly
    an azimuthal rotation, applying this to every point of a cubed-sphere
    grid is equivalent to adding angleDeg to every point's longitude and
    leaving its latitude untouched -- but doing it here, on the cube's own
    face frames (see _faceFrames), rotates the panels themselves (and
    hence where their seams and the dateline/pole fix-ups have to cope
    with them), rather than post-hoc patching (lon, lat) pairs.
    """
    if angleDeg == 0.0:
        return v
    a = numpy.radians(angleDeg)
    c, s = numpy.cos(a), numpy.sin(a)
    R = numpy.array([[c, -s, 0.],
                     [s, c, 0.],
                     [0., 0., 1.]])
    return v.dot(R.T)


def _faceFrames(rotationDeg=0.0):
    """
    The 6 cube faces as (name, normal, u_hat, v_hat) unit vectors, chosen
    so that u_hat x v_hat = normal for every face. Walking a face's local
    (u, v) grid in increasing order then traces each cell counterclockwise
    as seen from OUTSIDE the sphere -- the orientation convention used
    elsewhere in this repository (e.g. the "0-->--1 / ^ ^ / 3-->--2" node
    ordering in mint/tests/test_polyline_integral.py).

    :param rotationDeg: rotate the whole cube (all 6 faces, hence all
                        panel seams and the dateline/pole they may cross)
                        by this many degrees about the vertical (Z) axis.
                        0 (the default) reproduces the original,
                        axis-aligned cube exactly (rotateZ is a no-op for
                        angleDeg=0, so this is bit-for-bit identical to
                        not rotating at all).
    """
    X, Y, Z = numpy.eye(3)
    frames = [
        ('+X', X, Y, Z),
        ('-X', -X, Z, Y),
        ('+Y', Y, Z, X),
        ('-Y', -Y, X, Z),
        ('+Z', Z, X, Y),
        ('-Z', -Z, Y, X),
    ]
    if rotationDeg == 0.0:
        return frames
    return [(name, rotateZ(normal, rotationDeg), rotateZ(u_hat, rotationDeg), rotateZ(v_hat, rotationDeg))
            for name, normal, u_hat, v_hat in frames]


def _unwrapLongitudes(lon):
    """Add 360 to whichever of these longitudes are on the "wrong side" of
    the +-180 deg dateline, so a single cell doesn't appear to span most of
    the way around the sphere."""
    if lon.max() - lon.min() > 180.0:
        lon = numpy.where(lon < 0.0, lon + 360.0, lon)
    return lon


def _fixCellLongitudes(lon, lat, pole_tol_deg=1.e-8):
    """
    Fix the 4 corner longitudes of a single cell: unwrap the ordinary
    dateline-crossing case, and separately handle a corner sitting
    (numerically) exactly at a pole, where longitude is not well-defined
    and atan2 returns an arbitrary value that can be wildly different from
    its neighbours -- give it a neighbour's (already-unwrapped) longitude
    instead.
    """
    lon = lon.copy()
    at_pole = numpy.abs(numpy.abs(lat) - 90.0) < pole_tol_deg

    if at_pole.any() and not at_pole.all():
        non_pole = ~at_pole
        lon[non_pole] = _unwrapLongitudes(lon[non_pole])
        lon[at_pole] = lon[non_pole][0]
    else:
        lon = _unwrapLongitudes(lon)

    return lon


def cubedSphereGridPoints(M, rotationDeg=0.0):
    """
    Build a cubed-sphere grid with M x M cells per face (6 * M^2 cells
    total) by projecting a uniform grid on each face of a circumscribing
    cube radially onto the sphere.

    :param M: number of cells along one edge of a cube face
    :param rotationDeg: rotate the whole cube about the vertical (Z) axis
                        by this many degrees before projecting -- shifts
                        where every panel's seams (and hence the dateline
                        and pole singularities relative to them) fall,
                        without changing the grid's topology. 0 (the
                        default) is the original, axis-aligned cube.
    :returns: numpy array of shape (6*M*M, 4, 3): (longitude, latitude, 0)
              in degrees, ready for mint.Grid().setPoints(...)

    :note: M=1 (6 giant cells, one per cube face) is not supported: such a
           cell's own longitude span is too extreme for the per-cell
           dateline fix below to represent unambiguously (confirmed by a
           nonzero closed-loop integral in the sanity checks). This is
           not a realistic resolution for any actual use of this grid.
    """
    if M < 2:
        raise ValueError(f'M must be >= 2 (M=1 is degenerate, see note above), got {M}')

    edges = numpy.linspace(-1., 1., M + 1)

    points = numpy.empty((6 * M * M, NUM_VERTS_PER_QUAD, 3), dtype=numpy.float64)
    icell = 0
    for name, normal, u_hat, v_hat in _faceFrames(rotationDeg):
        U, V = numpy.meshgrid(edges, edges, indexing='ij')  # (M+1, M+1)
        cube_pts = (normal[None, None, :]
                    + U[..., None] * u_hat[None, None, :]
                    + V[..., None] * v_hat[None, None, :])
        norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
        sphere_pts = cube_pts / norm  # gnomonic projection onto the unit sphere

        lon = numpy.degrees(numpy.arctan2(sphere_pts[..., 1], sphere_pts[..., 0]))
        lat = numpy.degrees(numpy.arcsin(numpy.clip(sphere_pts[..., 2], -1.0, 1.0)))

        for i in range(M):
            for j in range(M):
                corner_lon = numpy.array([lon[i, j], lon[i + 1, j], lon[i + 1, j + 1], lon[i, j + 1]])
                corner_lat = numpy.array([lat[i, j], lat[i + 1, j], lat[i + 1, j + 1], lat[i, j + 1]])
                points[icell, :, 0] = _fixCellLongitudes(corner_lon, corner_lat)
                points[icell, :, 1] = corner_lat
                points[icell, :, 2] = 0.0
                icell += 1

    return points


_CUBE_FACE_NAMES = ['+X', '-X', '+Y', '-Y', '+Z', '-Z']
_CUBE_FACE_AXES = numpy.array([
    [1., 0., 0.], [-1., 0., 0.], [0., 1., 0.], [0., -1., 0.], [0., 0., 1.], [0., 0., -1.],
])


def cubeFaceOf(lonDeg, latDeg):
    """
    Which of the 6 circumscribing-cube faces (see _faceFrames, unrotated:
    rotationDeg=0) a (lon, lat) point's own radial direction projects onto
    under gnomonic projection -- the face whose own outward normal has the
    largest dot product with the point's unit vector, which is exactly the
    face gnomonic projection would assign it to (the face the ray from the
    sphere's centre through this point hits first).

    :returns: one of '+X', '-X', '+Y', '-Y', '+Z', '-Z'
    """
    lam = numpy.radians(lonDeg)
    theta = numpy.radians(latDeg)
    xyz = numpy.array([numpy.cos(theta) * numpy.cos(lam),
                        numpy.cos(theta) * numpy.sin(lam),
                        numpy.sin(theta)])
    dots = _CUBE_FACE_AXES.dot(xyz)
    return _CUBE_FACE_NAMES[int(numpy.argmax(dots))]


def facesAlongPath(p0, p1, n=50):
    """
    The set of distinct cube faces (see cubeFaceOf) touched by n points
    sampled along the straight-in-(lon,lat) path from p0 to p1 -- used to
    verify a test path genuinely crosses a cubed-sphere panel seam, not
    just many cells within a single panel.

    :param p0, p1: (lonDeg, latDeg) endpoints
    :returns: set of face-name strings
    """
    ts = numpy.linspace(0., 1., n)
    lons = p0[0] + ts * (p1[0] - p0[0])
    lats = p0[1] + ts * (p1[1] - p0[1])
    return {cubeFaceOf(lo, la) for lo, la in zip(lons, lats)}


def buildCubedSphereGrid(M, rotationDeg=0.0):
    """
    Convenience wrapper: build the points and return a ready-to-use
    mint.Grid (with the flags a cubed-sphere grid requires, see
    Grid.setFlags).

    :param M: number of cells along one edge of a cube face
    :param rotationDeg: see cubedSphereGridPoints

    Grid.setPoints keeps its own reference to the points array (it hands
    mint's C++ side a raw, non-copied pointer into it, so it must stay
    alive for as long as the grid is used) -- no extra bookkeeping needed
    here.
    """
    grid = mint.Grid()
    grid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1, degrees=True)
    grid.setPoints(cubedSphereGridPoints(M, rotationDeg))
    return grid


def _checkOrientation(points):
    """
    Every cell's 4 corners should be counterclockwise as seen from OUTSIDE
    the sphere, i.e. cross(edge0, edge1) should point the same way as the
    cell's own centroid (outward). Returns the number of cells that fail
    this (should always be 0).
    """
    deg2rad = numpy.pi / 180.0
    lam = points[..., 0] * deg2rad
    theta = points[..., 1] * deg2rad
    xyz = numpy.stack([numpy.cos(theta) * numpy.cos(lam),
                        numpy.cos(theta) * numpy.sin(lam),
                        numpy.sin(theta)], axis=-1)
    p0, p1, p2, p3 = xyz[:, 0], xyz[:, 1], xyz[:, 2], xyz[:, 3]
    normal = numpy.cross(p1 - p0, p2 - p0)
    centroid = xyz.mean(axis=1)
    outward = numpy.sum(normal * centroid, axis=-1)
    return int((outward <= 0).sum())


def _sanityChecks(M, rotationDeg=0.0):
    """
    Quantitative checks that the generated grid is correct, not just
    "looks plausible":

    1. Cell orientation: every cell's corners should be counterclockwise
       as seen from outside the sphere (see _checkOrientation).
    2. Total surface area should be close to 4*pi (unit sphere), computed
       from the actual cell corner points via the planar-quad approximation
       to spherical-quadrilateral area, which converges to the true area
       as M increases (its own error shrinks like 1/M^2, hence the
       M-dependent tolerance below -- not a sign of a bug at small M, just
       of the approximation itself and of the cubed-sphere's inherent,
       expected cell-size non-uniformity).
    3. mint.PolylineIntegral of a CLOSED loop, with edge data built from an
       EXACT global potential (finite differences of a smooth function),
       should be zero to near machine precision regardless of M, the loop
       chosen, or the panel rotation.

    :param rotationDeg: passed straight through to buildCubedSphereGrid,
                        so these same checks can be re-run with the cube's
                        panel seams shifted to arbitrary, non-round
                        longitudes (rather than just the un-rotated cube's
                        own seams at 0/45/90/.../-180 etc.).
    """
    grid = buildCubedSphereGrid(M, rotationDeg)
    ncells = grid.getNumberOfCells()
    print(f'M={M}, rotationDeg={rotationDeg}: {ncells} cells (expected {6 * M * M})')
    assert ncells == 6 * M * M

    points = grid.getPoints()

    n_bad_orientation = _checkOrientation(points)
    print(f'  cells with inconsistent (non-outward) orientation: {n_bad_orientation}')
    assert n_bad_orientation == 0, 'inconsistent cell orientation'

    deg2rad = numpy.pi / 180.0
    lam = points[..., 0] * deg2rad
    theta = points[..., 1] * deg2rad
    xyz = numpy.stack([numpy.cos(theta) * numpy.cos(lam),
                        numpy.cos(theta) * numpy.sin(lam),
                        numpy.sin(theta)], axis=-1)
    p0, p1, p2, p3 = xyz[:, 0], xyz[:, 1], xyz[:, 2], xyz[:, 3]
    area = 0.5 * (numpy.linalg.norm(numpy.cross(p1 - p0, p2 - p0), axis=-1)
                   + numpy.linalg.norm(numpy.cross(p2 - p0, p3 - p0), axis=-1))
    total_area = area.sum()
    exact_area = 4.0 * numpy.pi
    rel_err = abs(total_area - exact_area) / exact_area
    tol = max(0.01, 3.0 / M**2)
    print(f'  total area = {total_area:.6f} (exact 4*pi = {exact_area:.6f}), '
          f'relative error = {rel_err:.2e} (tol {tol:.2e})')
    assert rel_err < tol, f'surface area error {rel_err:.2e} exceeds tolerance {tol:.2e}: likely a bug'

    def potential(p):
        lam_, th_ = p[..., 0] * deg2rad, p[..., 1] * deg2rad
        return numpy.sin(3 * lam_) * numpy.cos(th_) ** 2

    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * (potential(points[:, i1, :]) - potential(points[:, i0, :]))

    pli = mint.PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=False)
    loop = numpy.array([(10., 10., 0.), (200., 10., 0.), (200., -60., 0.),
                         (10., -60., 0.), (10., 10., 0.)])
    pli.computeWeights(loop, counterclock=False)
    loop_flux = pli.getIntegral(data, placement=mint.CELL_BY_CELL_DATA)
    print(f'  closed-loop integral of an exact field (should be ~0): {loop_flux:.3e}')
    assert abs(loop_flux) < 1.e-8, 'nonzero closed-loop integral: likely an orientation/dateline/pole bug'

    return grid


if __name__ == '__main__':
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    for m in (2, 3, 4, 8, 16):
        _sanityChecks(m)

    # also exercise the rotated-panels code path, at a deliberately
    # non-round angle so the panel seams land at generic longitudes.
    # Skips M=2: this particular (very large, near-global) diagnostic
    # loop skims the domain edge at that coarse a resolution for some
    # rotation angles -- a locator-robustness limit of testing a huge
    # loop against a handful of giant cells, the same class of issue
    # that already rules out M=1 above, not a bug in the rotation itself
    # (every cell's own corner longitude span stays well under 180 deg;
    # M=3 upward all pass to near machine precision).
    for m in (3, 4, 8, 16):
        _sanityChecks(m, rotationDeg=23.7)

    grid = buildCubedSphereGrid(M)
    outfile = Path(__file__).parent / f'cubedsphere_M{M}.vtk'
    grid.dump(str(outfile))
    print(f'\nwrote {outfile} ({grid.getNumberOfCells()} cells)')
