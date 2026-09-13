#!/usr/bin/env python
"""
Validation for the getTangentVectors / isPointInQuad fix (see mntVectorInterp.h,
vmtCellLocator.cpp): build the SAME pole-containing cubed-sphere panel ('+Z')
that showed degraded (~O(h^1.9)) convergence when stored as (lon, lat, elev),
but this time store its points as genuine embedded (x, y, z) unit-sphere
Cartesian coordinates -- no lon/lat at all, so there is no coordinate
singularity at the pole in this representation -- and check whether
mint.VectorInterp's REAL (C++, just-rebuilt) reconstruction now converges at
the clean, non-degraded rate instead.

Field: same pole-asymmetric vector potential u = curl(psi * r_hat),
psi(lambda, theta) = cos(theta)*(1+sin(theta))*cos(lambda), reused from
pole_asymmetric_convergence.py, but evaluated here as a genuine AMBIENT (x,y,z)
vector: since p(lambda, theta) = (cos(theta)cos(lambda), cos(theta)sin(lambda),
sin(theta)) parametrizes the unit sphere, dp/dlambda and dp/dtheta are the
ambient-space tangent vectors (NOT unit -- |dp/dlambda| = cos(theta)), and
v_xyz = u_lam * dp/dlambda + u_theta * dp/dtheta is the true physical vector in
ambient (x,y,z) components, no metric/DEG2RAD correction needed (unlike the
(lon,lat)-basis case) -- this IS the raw coordinate-basis output an
xyz-embedded grid's VectorInterp should reproduce directly.

Edge circulations and the "exact" pointwise-reconstruction target are built
with whitney2.py's plain Gauss-Legendre quadrature helpers (coordinate-agnostic
-- they already work on any 2- or 3-component point representation), reusing
the SAME machinery as the earlier k=2 Whitney convergence scripts, but this
script only exercises the actual C++ mint.VectorInterp (k=1) reconstruction --
not the from-scratch Python one in whitney2.py.

Usage: python scripts/xyz_pole_panel_validation.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import mint

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from whitney2 import edgeIntegral  # noqa: E402
from pole_asymmetric_convergence import uDotDl  # noqa: E402

NUM_VERTS_PER_QUAD = 4


# ---------------------------------------------------------------------------
# Genuine ambient (x, y, z) vector field: v_xyz(p) = u_lam * dp/dlambda +
# u_theta * dp/dtheta, both tangent vectors evaluated at p's own (lambda, theta).
# ---------------------------------------------------------------------------
def vectorFieldXYZ(pts):
    """(..., 3) unit-sphere xyz points -> (..., 3) ambient physical vector."""
    x, y, z = pts[..., 0], pts[..., 1], pts[..., 2]
    lam = numpy.arctan2(y, x)
    the = numpy.arcsin(numpy.clip(z, -1.0, 1.0))

    u_lam = uDotDl(lam, the, 1.0, 0.0)
    u_the = uDotDl(lam, the, 0.0, 1.0)

    dpdlam = numpy.stack([-numpy.cos(the) * numpy.sin(lam),
                          numpy.cos(the) * numpy.cos(lam),
                          numpy.zeros_like(lam)], axis=-1)
    dpdthe = numpy.stack([-numpy.sin(the) * numpy.cos(lam),
                          -numpy.sin(the) * numpy.sin(lam),
                          numpy.cos(the)], axis=-1)

    return u_lam[..., None] * dpdlam + u_the[..., None] * dpdthe


# ---------------------------------------------------------------------------
# '+Z' panel only, stored as genuine xyz (no lon/lat conversion at all) --
# adapted from generate_cubedsphere_grid.cubedSphereGridPoints, keeping the
# unit-sphere xyz points it already computes internally instead of converting
# them to (lon, lat).
# ---------------------------------------------------------------------------
def buildPlusZPanelXYZ(M):
    """Returns (M*M, 4, 3) corners, genuine unit-sphere (x, y, z)."""
    normal, u_hat, v_hat = numpy.array([0., 0., 1.]), numpy.array([1., 0., 0.]), numpy.array([0., 1., 0.])
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = numpy.meshgrid(edges, edges, indexing='ij')
    cube_pts = (normal[None, None, :] + U[..., None] * u_hat[None, None, :]
                + V[..., None] * v_hat[None, None, :])
    norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
    sphere_pts = cube_pts / norm  # (M+1, M+1, 3) unit-sphere xyz -- the pole is at U=V=0, i.e. cell centre

    points = numpy.empty((M * M, NUM_VERTS_PER_QUAD, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            points[icell, 0] = sphere_pts[i, j]
            points[icell, 1] = sphere_pts[i + 1, j]
            points[icell, 2] = sphere_pts[i + 1, j + 1]
            points[icell, 3] = sphere_pts[i, j + 1]
            icell += 1
    return points


def buildExactEdgeData(points, n=16):
    """Cell-by-cell circulations, in raw xyz coordinates (no lon/lat)."""
    ncells = points.shape[0]
    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        for icell in range(ncells):
            data[icell, i0] = sign * edgeIntegral(points[icell, i0], points[icell, i1], vectorFieldXYZ, None, n)
    return data


def buildXYZGrid(points):
    """mint.Grid from genuine xyz points -- NOT lon/lat: fixLonAcrossDateline,
    averageLonAtPole and degrees are all disabled (there is no dateline/pole
    special-casing to apply to an ambient xyz representation)."""
    grid = mint.Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=0)
    grid.setPoints(points)
    return grid


XI_ETA_POINTS = [
    ('cell centre (0.50, 0.50)', 0.50, 0.50),
    ('off-centre  (0.25, 0.75)', 0.25, 0.75),
]
RESOLUTIONS = (4, 8, 16, 32, 64)


def bilinearCellTargets(points, xi, eta):
    p0, p1, p2, p3 = points[:, 0], points[:, 1], points[:, 2], points[:, 3]
    return ((1 - xi) * (1 - eta) * p0 + xi * (1 - eta) * p1
             + xi * eta * p2 + (1 - xi) * eta * p3)


def computeErrors(M):
    points = buildPlusZPanelXYZ(M)
    data = buildExactEdgeData(points)

    errors = {}
    for label, xi, eta in XI_ETA_POINTS:
        targets = bilinearCellTargets(points, xi, eta)
        # renormalize onto the unit sphere -- bilinear interpolation of 4
        # unit vectors is not itself unit length
        targets = targets / numpy.linalg.norm(targets, axis=-1, keepdims=True)

        grid = buildXYZGrid(points)
        vi = mint.VectorInterp()
        vi.setGrid(grid)
        vi.buildLocator(numCellsPerBucket=128, periodX=0.0, enableFolding=False)
        numBad = vi.findPoints(targets, tol2=1.e-6)
        assert numBad == 0, f'M={M}: {numBad}/{len(targets)} target points at {label} not found'

        numeric = vi.getEdgeVectors(data, placement=mint.CELL_BY_CELL_DATA)
        exact = vectorFieldXYZ(targets)
        errors[label] = numpy.linalg.norm(numeric - exact, axis=-1)
    return errors


def fitAlpha(Ms, errs):
    Ms = numpy.asarray(Ms, dtype=numpy.float64)
    errs = numpy.asarray(errs, dtype=numpy.float64)
    slope, intercept = numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)
    return -slope, numpy.exp(intercept)


def main():
    rmsByLabel = {label: [] for label, _, _ in XI_ETA_POINTS}
    for M in RESOLUTIONS:
        errors = computeErrors(M)
        pieces = []
        for label, _, _ in XI_ETA_POINTS:
            rms = float(numpy.sqrt((errors[label] ** 2).mean()))
            rmsByLabel[label].append(rms)
            pieces.append(f'{label.split("(")[0].strip()}={rms:.3e}')
        print(f'M={M:3d} (ncells={M * M:5d}): ' + ', '.join(pieces))

    fig, ax = plt.subplots(figsize=(7, 6))
    Ms = numpy.asarray(RESOLUTIONS, dtype=numpy.float64)
    for label, _, _ in XI_ETA_POINTS:
        alpha, prefactor = fitAlpha(RESOLUTIONS, rmsByLabel[label])
        ax.loglog(Ms, rmsByLabel[label], 'o-', label=f'{label} (alpha={alpha:.2f})')
        print(f'{label}: fitted error ~ M^-{alpha:.3f}')

    ax.set_xlabel('M (cells per panel edge)')
    ax.set_ylabel('RMS |v_numeric - v_exact| (ambient xyz)')
    ax.set_title("mint.VectorInterp on the '+Z' (pole) panel, genuine xyz points\n"
                  "(post-fix: no lon/lat coordinate singularity at the pole)")
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    outfile = HERE / 'xyz_pole_panel_validation.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
