#!/usr/bin/env python
"""
Python prototype + validation for extending mint's RegridEdges to genuine
(x, y, z) embedded data (see the conversation this implements: RegridEdges'
line-crossing machinery, PolysegmentIter -> vmtLonLatCellLocator::
findIntersectionsWithLine -> collectIntersectionPoints, is 2D-only --
LineLineIntersector solves an exact line-vs-line intersection, which is
well-posed when source/destination edges share a flat (lon, lat) plane but
NOT for a genuinely curved 3D surface).

Key finding this script starts from (verified numerically before writing
any of this): treating a destination edge as a straight 3D chord (as
instructed) does NOT mean "cast a ray and see which cell surfaces it
pierces" -- a chord between two points on a curved surface (e.g. a sphere)
dips BELOW the surface everywhere except at its own endpoints (e.g. ~0.0038
units below a unit sphere at the midpoint of a 10-degree-long chord), so a
literal ray-vs-surface intersection test would find zero intermediate
cells for any chord spanning more than one cell width -- exactly the
"coarse destination edge crossing several finer source cells" case
regridding actually needs to handle. So the algorithm here is NOT built on
VTK's IntersectWithLine (ray-cast); it's point-sampling + bisection:

  1. Locate the cell containing the chord's current point via the SAME
     Newton/Gauss-Newton bilinear inversion vmtLonLatCellLocator's
     invertSphericalBilinearPatch already uses for the spherical case
     (locateInCell below), just for a flat (non-spherical) bilinear patch.
  2. Bisect forward along the chord to find where the point exits that
     cell's parametric [0,1]^2 domain.
  3. Repeat from the exit point for the next cell, until the whole [0,1]
     line-parameter range is covered.

This produces exactly the (cellId, ta, tb) segment list PolysegmentIter
needs, from a totally different (but here, actually correct) geometric
primitive than the ray-cast approach.

Validated two ways below: a flat (z=0) Cartesian case (isolates whether the
segment-walking bookkeeping itself is correct, no curved-surface error to
confound it) and a genuinely curved xyz cubed-sphere case (also picks up
the "chord approximates the true edge" error, expected to shrink with
resolution the same way mint's existing lon-lat flat-chord approximation
does) -- both compared against the EXACT analytic edge circulation (not
against the C++ pipeline, since that's lon-lat-only today).

Usage: python scripts/xyz_regrid_prototype.py
"""
import sys
from pathlib import Path

import numpy

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from whitney2 import bilinearMap, bilinearPartials, edgeIntegral, gaussLegendre01  # noqa: E402


# ---------------------------------------------------------------------------
# 1. Locate a point within a single bilinear cell (Gauss-Newton, same recipe
#    as vmtLonLatCellLocator::invertSphericalBilinearPatch, but for a flat
#    bilinear patch -- no slerp/spherical assumption).
# ---------------------------------------------------------------------------
def locateInCell(corners, target, xi0=0.5, eta0=0.5, maxIter=30, tol2=1.e-24, h=1.e-6):
    """
    :param corners: (4, 2 or 3) cell corners
    :param target: (2 or 3,) point to locate
    :returns: (xi, eta, dist2) -- dist2 is the squared distance from the
              bilinear patch's closest point (to xi,eta) to target; 0 for a
              flat 2D case where target is exactly reachable, generally > 0
              for a 3D point near (not exactly on) a curved patch.
    """
    xi, eta = xi0, eta0
    for _ in range(maxIter):
        X = bilinearMap(corners, numpy.array([xi]), numpy.array([eta]))[0]
        res = X - target
        dXdxi, dXdeta = bilinearPartials(corners, numpy.array([xi]), numpy.array([eta]))
        dXdxi, dXdeta = dXdxi[0], dXdeta[0]
        a11 = numpy.dot(dXdxi, dXdxi)
        a12 = numpy.dot(dXdxi, dXdeta)
        a22 = numpy.dot(dXdeta, dXdeta)
        b1 = -numpy.dot(dXdxi, res)
        b2 = -numpy.dot(dXdeta, res)
        det = a11 * a22 - a12 * a12
        if abs(det) < 1.e-300:
            break
        dxi = (a22 * b1 - a12 * b2) / det
        deta = (a11 * b2 - a12 * b1) / det
        xi += dxi
        eta += deta
        if dxi * dxi + deta * deta < tol2:
            break
    X = bilinearMap(corners, numpy.array([xi]), numpy.array([eta]))[0]
    dist2 = float(numpy.dot(X - target, X - target))
    return xi, eta, dist2


def findContainingCell(allCorners, target, tol=1.e-6, xi0=0.5, eta0=0.5):
    """Linear scan over all cells (fine for a small prototype grid; a real
    implementation would use the bucket-indexed locator, exactly as
    vmtXYZCellLocator's FindCell already does for point queries)."""
    best = (None, None, None, numpy.inf)
    for icell in range(allCorners.shape[0]):
        xi, eta, dist2 = locateInCell(allCorners[icell], target, xi0, eta0)
        if -tol <= xi <= 1 + tol and -tol <= eta <= 1 + tol and dist2 < best[3]:
            best = (icell, xi, eta, dist2)
    return best  # (cellId or None, xi, eta, dist2)


# ---------------------------------------------------------------------------
# 2. Walk the chord p0->p1, bisecting at each cell exit -- the actual
#    xyz-native replacement for findIntersectionsWithLine.
# ---------------------------------------------------------------------------
def findIntersectionsWithLineXYZ(allCorners, p0, p1, tol=1.e-6, nBisect=60):
    """
    :returns: list of (cellId, ta, tb) covering [0, 1] (ta, tb the linear
              parameter along p0->p1 where this segment enters/exits cellId)
    """
    p0 = numpy.asarray(p0, dtype=numpy.float64)
    p1 = numpy.asarray(p1, dtype=numpy.float64)
    segments = []
    t = 0.0
    xiGuess, etaGuess = 0.5, 0.5
    while t < 1.0 - 1.e-12:
        p = p0 + t * (p1 - p0)
        cellId, xi, eta, dist2 = findContainingCell(allCorners, p, tol, xiGuess, etaGuess)
        if cellId is None:
            # chord left the grid entirely -- stop (a real implementation
            # would report this as "outside", matching FindCell's -1)
            break

        def stillInside(tt):
            pp = p0 + tt * (p1 - p0)
            xi_, eta_, _ = locateInCell(allCorners[cellId], pp, xi, eta)
            return -tol <= xi_ <= 1 + tol and -tol <= eta_ <= 1 + tol

        if stillInside(1.0):
            tb = 1.0
        else:
            lo, hi = t, 1.0
            for _ in range(nBisect):
                mid = 0.5 * (lo + hi)
                if stillInside(mid):
                    lo = mid
                else:
                    hi = mid
            tb = lo

        segments.append((cellId, t, tb))
        xiGuess, etaGuess = xi, eta
        t = tb + 1.e-9  # nudge past the boundary
    return segments


# ---------------------------------------------------------------------------
# 3. computeWeight, ported unchanged from mntWeights.cpp (already
#    coordinate-agnostic -- pure parametric-space math, verified by reading
#    the C++ source) -- reused as-is, not reinvented.
# ---------------------------------------------------------------------------
def computeWeight(srcXi0, srcXi1, xia, xib):
    weight = 1.0
    sgn = 0.0
    for d in range(2):
        xiM = 0.5 * (xia[d] + xib[d])
        dxi = xib[d] - xia[d]
        xm = 0.5 * (srcXi0[d] + srcXi1[d])
        sgn += srcXi1[d] - srcXi0[d]
        xm00, xm05, xm10 = xm, xm - 0.5, xm - 1.0
        lag00 = 2.0 * xm05 * xm10
        lag05 = -4.0 * xm00 * xm10
        lag10 = 2.0 * xm00 * xm05
        weight *= (1.0 - xiM) * lag00 + dxi * lag05 + xiM * lag10
    return sgn * weight


REF_EDGE_PARAM_COORDS = numpy.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])  # corner -> (xi, eta)


def regridEdgeXYZ(srcCorners, srcEdgeData, dstP0, dstP1, tol=1.e-6):
    """
    The actual regridding computation for ONE destination edge: find the
    source cells/segments the chord dstP0->dstP1 crosses, and accumulate
    each crossed source cell's own 4 edge circulations, weighted by
    computeWeight -- this is exactly RegridEdges' inner loop
    (mnt_regridedges_computeWeights), just driven by findIntersectionsWithLineXYZ
    instead of PolysegmentIter/findIntersectionsWithLine.

    :param srcCorners: (ncells, 4, 2 or 3) source grid corners
    :param srcEdgeData: (ncells, 4) source edge circulations (cell-by-cell)
    :param dstP0, dstP1: destination edge endpoints
    :returns: regridded circulation for this destination edge
    """
    segments = findIntersectionsWithLineXYZ(srcCorners, dstP0, dstP1, tol=tol)
    total = 0.0
    for cellId, ta, tb in segments:
        corners = srcCorners[cellId]
        pa = numpy.asarray(dstP0) + ta * (numpy.asarray(dstP1) - numpy.asarray(dstP0))
        pb = numpy.asarray(dstP0) + tb * (numpy.asarray(dstP1) - numpy.asarray(dstP0))
        xiA0, xiA1, _ = locateInCell(corners, pa)
        xiB0, xiB1, _ = locateInCell(corners, pb)
        xia = (xiA0, xiA1)
        xib = (xiB0, xiB1)
        for srcEdgeIndex in range(4):
            is0, is1 = srcEdgeIndex, (srcEdgeIndex + 1) % 4
            weight = computeWeight(REF_EDGE_PARAM_COORDS[is0], REF_EDGE_PARAM_COORDS[is1], xia, xib)
            total += weight * srcEdgeData[cellId, srcEdgeIndex]
    return total


# ---------------------------------------------------------------------------
# Validation 1: flat (z=0) Cartesian grid, destination COARSER than source
# (so a single destination edge spans several source cells -- the case that
# actually stresses the segment-walking bookkeeping, and that a literal
# ray-cast would get wrong). No curvature here, so this isolates whether
# the walking/bisection/weight-accumulation logic itself is correct from
# any curved-surface chord-approximation error (that's validation 2).
# ---------------------------------------------------------------------------
KX, KY = 2.0 * numpy.pi, numpy.pi


def vectorFieldFlat(pts):
    """u = dpsi/dy, w = -dpsi/dx for psi(x,y) = sin(KX*x)*sin(KY*y) (same
    field as whitney2_convergence_cartesian.py) -- genuinely rotational,
    not a gradient, so edge circulations have to be actually integrated."""
    x, y = pts[..., 0], pts[..., 1]
    u = KY * numpy.sin(KX * x) * numpy.cos(KY * y)
    w = -KX * numpy.cos(KX * x) * numpy.sin(KY * y)
    return numpy.stack([u, w, numpy.zeros_like(u)], axis=-1)


def buildCartesianGrid(M, Lx=1.0, Ly=1.0):
    xs = numpy.linspace(0., Lx, M + 1)
    ys = numpy.linspace(0., Ly, M + 1)
    corners = numpy.empty((M * M, 4, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            corners[icell, 0] = (xs[i], ys[j], 0.)
            corners[icell, 1] = (xs[i + 1], ys[j], 0.)
            corners[icell, 2] = (xs[i + 1], ys[j + 1], 0.)
            corners[icell, 3] = (xs[i], ys[j + 1], 0.)
            icell += 1
    return corners


def buildExactEdgeData(corners, field, n=16):
    """Cell-by-cell raw (corner_i -> corner_{i+1}) circulations -- NOT the
    sign-flipped convention mnt_vectorinterp uses for its bilinear blend;
    RegridEdges reads vtkCell's own parametric corner order directly (see
    mntRegridEdges.cpp), so this is the convention to match here."""
    ncells = corners.shape[0]
    data = numpy.zeros((ncells, 4))
    for i0 in range(4):
        i1 = (i0 + 1) % 4
        for icell in range(ncells):
            data[icell, i0] = edgeIntegral(corners[icell, i0], corners[icell, i1], field, None, n)
    return data


def validateFlatCase(Msrc, Mdst):
    srcCorners = buildCartesianGrid(Msrc)
    dstCorners = buildCartesianGrid(Mdst)
    srcEdgeData = buildExactEdgeData(srcCorners, vectorFieldFlat)

    errs = []
    for icell in range(dstCorners.shape[0]):
        c = dstCorners[icell]
        for i0 in range(4):
            i1 = (i0 + 1) % 4
            exact = edgeIntegral(c[i0], c[i1], vectorFieldFlat, None, n=16)
            regridded = regridEdgeXYZ(srcCorners, srcEdgeData, c[i0], c[i1])
            errs.append(regridded - exact)
    errs = numpy.asarray(errs)
    return float(numpy.sqrt(numpy.mean(errs ** 2))), float(numpy.max(numpy.abs(errs)))


# ---------------------------------------------------------------------------
# Validation 2: genuinely curved xyz cubed-sphere panel, same idea, now also
# picking up the "chord approximates the true edge" error -- expected to
# behave like mint's existing lon-lat flat-chord approximation: fine for
# reasonable resolutions, degrading at some point for a fixed dst/src
# resolution ratio when the chord's sag stops being small relative to cell
# size. Per the user: that eventual breakdown is expected and acceptable,
# not something this prototype needs to work around.
# ---------------------------------------------------------------------------
def buildPlusXPanelXYZ(M):
    normal, u_hat, v_hat = numpy.array([1., 0., 0.]), numpy.array([0., 1., 0.]), numpy.array([0., 0., 1.])
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = numpy.meshgrid(edges, edges, indexing='ij')
    cube_pts = normal[None, None, :] + U[..., None] * u_hat[None, None, :] + V[..., None] * v_hat[None, None, :]
    norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
    sphere_pts = cube_pts / norm
    corners = numpy.empty((M * M, 4, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            corners[icell, 0] = sphere_pts[i, j]
            corners[icell, 1] = sphere_pts[i + 1, j]
            corners[icell, 2] = sphere_pts[i + 1, j + 1]
            corners[icell, 3] = sphere_pts[i, j + 1]
            icell += 1
    return corners


def vectorFieldXYZ(pts):
    """Same pole-asymmetric vector potential used throughout this codebase
    (pole_asymmetric_convergence.py), reused via
    xyz_pole_panel_validation.vectorFieldXYZ."""
    from xyz_pole_panel_validation import vectorFieldXYZ as _vf
    return _vf(pts)


def validateCurvedCase(Msrc, Mdst):
    srcCorners = buildPlusXPanelXYZ(Msrc)
    dstCorners = buildPlusXPanelXYZ(Mdst)
    srcEdgeData = buildExactEdgeData(srcCorners, vectorFieldXYZ)

    errs = []
    for icell in range(dstCorners.shape[0]):
        c = dstCorners[icell]
        for i0 in range(4):
            i1 = (i0 + 1) % 4
            exact = edgeIntegral(c[i0], c[i1], vectorFieldXYZ, None, n=16)
            regridded = regridEdgeXYZ(srcCorners, srcEdgeData, c[i0], c[i1])
            errs.append(regridded - exact)
    errs = numpy.asarray(errs)
    return float(numpy.sqrt(numpy.mean(errs ** 2))), float(numpy.max(numpy.abs(errs)))


if __name__ == '__main__':
    print('=== flat Cartesian, destination coarser than source ===')
    for Msrc, Mdst in [(4, 2), (8, 2), (8, 4), (16, 4)]:
        rms, mx = validateFlatCase(Msrc, Mdst)
        print(f'  src={Msrc:2d}x{Msrc:<2d} dst={Mdst:2d}x{Mdst:<2d}: rms={rms:.3e}  max={mx:.3e}')

    print('\n=== curved xyz cubed-sphere panel (+X), destination coarser than source ===')
    for Msrc, Mdst in [(4, 2), (8, 2), (8, 4), (16, 4)]:
        rms, mx = validateCurvedCase(Msrc, Mdst)
        print(f'  src={Msrc:2d}x{Msrc:<2d} dst={Mdst:2d}x{Mdst:<2d}: rms={rms:.3e}  max={mx:.3e}')
