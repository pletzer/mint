#!/usr/bin/env python
"""
Convergence of mint.VectorInterp's vector reconstruction (from edge-integral,
1-form data) on cubed-sphere grids of increasing resolution M, for target
points placed at arbitrary parametric locations (xi, eta) inside each cell --
not just cell centres.

The vector field is derived from a VECTOR potential, not a scalar one: as in
work/vector_potential/formulas.py ("Velocity from a vector potential on the
sphere"),

    u = curl(psi * r_hat),   psi(lambda, theta) = cos(theta)*(1+sin(theta))*cos(lambda)

i.e. psi*r_hat (a purely radial vector field) is the vector potential, and u
is its curl -- the same pole-asymmetric wavenumber-1 stream function already
used and sympy-verified in scripts/pole_asymmetric_convergence.py (finite
nonzero velocity at the north pole, exactly zero at the south pole). This is
deliberately NOT a gradient/scalar-potential field like the one in
mint/tests/test_vector_interp_cubedsphere.py: a gradient field's edge
integrals are exact by the fundamental theorem of calculus (just endpoint
differences), which is too special a case to stress the reconstruction at
generic interior points. A genuine (non-exact) 1-form has real circulation,
so its edge data has to come from actually integrating u.dl along each edge.

That edge integral is computed EXACTLY here (buildExactEdgeData below, using
pole_asymmetric_convergence.exactLineIntegral per edge, sympy, no
quadrature) rather than via that same module's straightLineIntegral (an
n-point numerical quadrature). Precisely what "exact" removes: both
functions integrate along the SAME path -- straight in (lon, lat) from
corner to corner -- exactLineIntegral just evaluates that integral in closed
form instead of approximating it, so this eliminates quadrature error only.
It does NOT touch the other, independent source of error already present in
either version: modelling a cubed-sphere edge as a flat (lon, lat) chord in
the first place, when the true mesh edge is closer to a great-circle arc.
Quadrature error was separately verified elsewhere (see
straightLineIntegral's and exactLineIntegral's own docstrings) to be ~1e-6
at the n=50 default, i.e. already negligible next to the errors under study
-- so switching to the exact version is expected to leave the fitted alpha
essentially unchanged; if it moved a lot instead, that would mean the
quadrature error was quietly contaminating the result after all. Both
remaining error sources -- the flat-chord edge-geometry approximation, and
the Whitney-1-form reconstruction itself -- are still lumped together in
what's measured below; separating those two would require exactly
integrating along the true great-circle edge, which is a much harder
(non-elementary, since the arc isn't linear in any angle parametrisation)
symbolic integral than exactLineIntegral solves, and is not attempted here.

exactLineIntegral redoes a fresh sympy symbolic integration for every single
edge, not vectorizable (see its own docstring for why); its docstring
estimates ~2-3ms/call, but measured directly for THIS field it's closer to
~60ms/call (evidently field-dependent -- more trig terms to expand/simplify
than whatever it was originally timed on). At that rate RESOLUTIONS below is
kept to (4, 8, 16) cells/panel-edge, not 32 or 64 (M=32 alone would be
6*32^2*4 ~= 24600 edge integrals, ~26 minutes on its own) -- enough points
for a fit and a direct, apples-to-apples comparison against the quadrature
version's M=4/8/16 numbers, at a still-substantial but practical ~8-9
minutes total runtime.

Parametric coordinates (xi, eta) in [0, 1]^2, per cell:

    3-->--2
    |     |          target(xi, eta) = (1-xi)(1-eta)*p0 + xi(1-eta)*p1
    ^     ^                          +      xi*eta *p2 + (1-xi)*eta *p3
    |     |
    0-->--1

i.e. plain bilinear interpolation of the cell's own 4 corners in (lon, lat)
space -- the same convention test_vector_interp_cubedsphere.py's
_cellCenters uses for (xi, eta) = (0.5, 0.5) (the average of the 4 corners),
generalised to arbitrary (xi, eta). This is a simple, self-consistent way to
place a controllable target point inside each cell; it is NOT necessarily
identical to mint's own internal parametric coordinate for a cubed-sphere
cell (a spherical bilinear patch / double-slerp, see vmtCellLocator -- the
pole-locator fix), so no assumption is made that this script's xi/eta lines
up exactly with what mint.VectorInterp would report as a target point's own
pcoords.

xi and eta must be STRICTLY inside (0, 1) (enforced by bilinearCellTargets),
for two separate reasons: (1) xi or eta = 0 or 1 puts the target exactly on
a shared cell edge (a shared corner if both are), where mint's FindCell
search can legitimately resolve to any of the neighbouring cells that also
touch that edge/corner -- which one it picks decides whose edge data the
vector is reconstructed from, so the "result" there is an implementation-
defined tie-break, not a meaningful sample of this scheme's error; (2) even
strictly inside a cell, this script's planar bilinear-in-(lon,lat) target is
only an approximation of mint's own (spherical-patch) cell geometry, so
xi/eta also need enough margin (see XI_ETA_POINTS) that this approximation
reliably still falls inside mint's own notion of the same cell.

For each (xi, eta) and each M, the RMS error over all 6*M^2 cells is fit to
error ~ M^{-alpha} via a log-log least-squares line; alpha is reported per
(xi, eta) location, plus one pooled over all sampled locations and cells
together.

Usage: python scripts/vector_potential_reconstruction_convergence.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import VectorInterp, CELL_BY_CELL_DATA, NUM_EDGES_PER_QUAD

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from generate_cubedsphere_grid import buildCubedSphereGrid  # noqa: E402
from pole_asymmetric_convergence import uDotDl, exactLineIntegral  # noqa: E402

DEG2RAD = numpy.pi / 180.0


# ---------------------------------------------------------------------------
# The exact vector field u = curl(psi * r_hat), in the same raw (dlambda,
# dtheta) coordinate-basis convention mint.VectorInterp reconstructs in (see
# test_vector_interp_cubedsphere.py's exactGradient for the analogous
# scalar-potential case): uDotDl(lam, theta, dlam, dtheta) is linear in
# (dlam, dtheta), so its two coefficients -- obtained by evaluating it at
# (dlam, dtheta) = (1, 0) and (0, 1) -- are exactly u's two coordinate-basis
# components.
# ---------------------------------------------------------------------------
def exactVector(lonDeg, latDeg):
    """
    Exact u = curl(psi * r_hat) at (lonDeg, latDeg), expressed in the grid's
    own (lon_deg, lat_deg) coordinate basis (matching mint.VectorInterp's raw
    output -- no spherical metric correction).

    uDotDl(lam, theta, dlam, dtheta) gives u_lam*dlam + u_theta*dtheta with
    lam/theta in RADIANS (needed for correct trig) and (u_lam, u_theta) the
    1-form's RADIAN-basis components, i.e. the coefficients of d(lam_rad),
    d(theta_rad). The edge-integral VALUES built by buildEdgeData below are
    already correct as is: a line integral is coordinate-basis-independent,
    and straightLineIntegral correctly uses the true radian angular
    displacement along each edge. But mint.VectorInterp's RECONSTRUCTED
    VECTOR is expressed w.r.t. the grid's raw, DEGREE-valued coordinate
    basis (it only ever sees the stored (lon_deg, lat_deg) numbers) -- so
    the target this is compared against needs the same chain-rule DEG2RAD
    factor as test_vector_interp_cubedsphere.py's exactGradient:
    d(lam_rad) = DEG2RAD * d(lon_deg) => degree-basis components are
    DEG2RAD * (radian-basis components), same physical 1-form, different
    basis. Omitting this factor was tried first and reproduced almost
    exactly |exact| as the error at every M (mint's output ends up ~57x
    smaller than an unscaled "exact", i.e. negligible next to it) --
    confirming this is the right fix, not a coincidence.

    Accepts scalar or array lonDeg/latDeg; returns an array of shape
    (..., 3) with a trailing zero z-component.
    """
    lam = lonDeg * DEG2RAD
    theta = latDeg * DEG2RAD
    u_lam = uDotDl(lam, theta, 1.0, 0.0)
    u_theta = uDotDl(lam, theta, 0.0, 1.0)
    return DEG2RAD * numpy.stack([u_lam, u_theta, numpy.zeros_like(u_lam)], axis=-1)


EDGE_MARGIN = 1.e-6  # keep target points strictly off cell edges/corners, see bilinearCellTargets


def bilinearCellTargets(points, xi, eta):
    """
    One target point per cell, at parametric (xi, eta) -- see module
    docstring for the corner ordering/convention.

    :param points: (ncells, 4, 3) cell corner points, as returned by
                    Grid.getPoints() (lon, lat, elevation, in degrees)
    :param xi, eta: scalars, STRICTLY inside (0, 1) (enforced below)
    :returns: (ncells, 2) array of (lon, lat) target points
    :raises ValueError: if xi or eta is 0, 1, or within EDGE_MARGIN of
                either -- such a point sits exactly on a shared cell
                edge/corner (a corner if BOTH are), which mint's FindCell
                search can legitimately resolve to any of the neighbouring
                cells that also touch that edge/corner. Which cell it picks
                then decides whose 4 edge-data values the vector gets
                reconstructed from, making the result an implementation-
                defined tie-break rather than a meaningful (xi, eta) sample
                -- not something a convergence-vs-M study should be
                sensitive to.
    """
    if not (EDGE_MARGIN < xi < 1. - EDGE_MARGIN) or not (EDGE_MARGIN < eta < 1. - EDGE_MARGIN):
        raise ValueError(
            f'(xi, eta) = ({xi}, {eta}) touches a cell boundary (both must be '
            f'strictly inside (0, 1)): a target point there lies exactly on a '
            f'shared edge/corner, where FindCell can resolve to any of several '
            f'neighbouring cells -- see bilinearCellTargets\' docstring.')
    p0, p1, p2, p3 = points[:, 0, :2], points[:, 1, :2], points[:, 2, :2], points[:, 3, :2]
    return ((1. - xi) * (1. - eta) * p0 + xi * (1. - eta) * p1
             + xi * eta * p2 + (1. - xi) * eta * p3)


def buildExactEdgeData(grid):
    """
    Cell-by-cell edge-integrated data for u = curl(psi * r_hat) on the
    grid's own edges, EXACT (sympy, no quadrature) -- see module docstring
    for exactly what "exact" does and doesn't remove here. Otherwise
    mirrors pole_asymmetric_convergence.buildEdgeData's structure (same
    corner/edge ordering, same sign convention for edges 2 and 3 running
    backwards relative to the (i0 -> i1) node order):

        3-->--2
        |     |
        ^     ^
        |     |
        0-->--1
    """
    ncells = grid.getNumberOfCells()
    points = grid.getPoints()  # (ncells, 4, 3): lon, lat, elev in degrees
    data = numpy.zeros((ncells, NUM_EDGES_PER_QUAD))
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        for icell in range(ncells):
            data[icell, i0] = sign * exactLineIntegral(points[icell, i0, :2], points[icell, i1, :2])
    return data


def computeErrors(M, xietaList):
    """
    Build the cubed-sphere grid and its exact edge data once for this M,
    then for each (label, xi, eta) in xietaList, locate the corresponding
    target point inside every cell and reconstruct the vector there.

    A fresh mint.VectorInterp is built (setGrid + buildLocator + findPoints)
    per (xi, eta), matching the pattern used everywhere else in this
    codebase (e.g. test_vector_interp_cubedsphere.py) rather than reusing
    one locator across multiple findPoints calls, which is not a pattern
    exercised/verified elsewhere.

    :returns: {label: err} where err is the (ncells,) array of
              |v_numeric - v_exact| for that (xi, eta), one value per cell
    """
    grid = buildCubedSphereGrid(M)
    data = buildExactEdgeData(grid)
    points = grid.getPoints()

    errors = {}
    for label, xi, eta in xietaList:
        targetsLonLat = bilinearCellTargets(points, xi, eta)
        targets = numpy.zeros((targetsLonLat.shape[0], 3))
        targets[:, :2] = targetsLonLat

        vi = VectorInterp()
        vi.setGrid(grid)
        vi.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
        numBad = vi.findPoints(targets, tol2=1.e-8)
        assert numBad == 0, (f'M={M}: {numBad} target points at {label} '
                              f'fell outside their own cell')

        numeric = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
        exact = exactVector(targetsLonLat[:, 0], targetsLonLat[:, 1])
        errors[label] = numpy.linalg.norm(numeric - exact, axis=-1)
    return errors


# ---------------------------------------------------------------------------
# Resolutions and the set of (xi, eta) target locations to probe. xi/eta are
# kept in [0.1, 0.9] -- see module docstring for why extreme values near a
# cell's own boundary are avoided.
# ---------------------------------------------------------------------------
RESOLUTIONS = (4, 8, 16)  # not 32/64: exactLineIntegral is per-edge sympy, see module docstring

XI_ETA_POINTS = [
    ('cell centre (0.50, 0.50)', 0.50, 0.50),
    ('off-centre  (0.25, 0.75)', 0.25, 0.75),
    ('near corner (0.10, 0.10)', 0.10, 0.10),
    ('near edge   (0.50, 0.10)', 0.50, 0.10),
]
POOLED_LABEL = 'pooled (all locations, all cells)'


def fitAlpha(Ms, errs):
    """Least-squares fit of log(err) = -alpha*log(M) + c; returns (alpha, prefactor)."""
    Ms = numpy.asarray(Ms, dtype=numpy.float64)
    errs = numpy.asarray(errs, dtype=numpy.float64)
    slope, intercept = numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)
    return -slope, numpy.exp(intercept)


def plotAndFitAlpha(rmsByLabel):
    """
    Log-log plot of RMS reconstruction error vs M, one curve per (xi, eta)
    location plus the pooled curve, each with its own power-law fit line.

    :param rmsByLabel: {label: (list of M, list of RMS error)}
    :returns: {label: alpha}
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.tab10.colors
    alphas = {}

    for i, (label, (Ms, errs)) in enumerate(rmsByLabel.items()):
        alpha, prefactor = fitAlpha(Ms, errs)
        alphas[label] = alpha

        color = colors[i % len(colors)]
        style = 'o-' if label != POOLED_LABEL else 's-'
        lw = 1.5 if label != POOLED_LABEL else 2.5
        ax.loglog(Ms, errs, style, color=color, linewidth=lw,
                   label=f'{label}  (alpha={alpha:.2f})')
        Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
        ax.loglog(Ms_arr, prefactor * Ms_arr ** (-alpha), '--', color=color, linewidth=1, alpha=0.5)

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel(r'RMS $|\mathbf{v}_{numeric} - \mathbf{v}_{exact}|$')
    ax.set_title('Vector reconstruction error vs cubed-sphere resolution\n'
                  r'$\mathbf{u} = \nabla\times(\psi\,\hat{r})$, '
                  r'$\psi=\cos\theta(1+\sin\theta)\cos\lambda$')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'vector_potential_reconstruction_convergence.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')
    return alphas


def main():
    rmsByLabel = {label: ([], []) for label, _, _ in XI_ETA_POINTS}
    rmsByLabel[POOLED_LABEL] = ([], [])

    for M in RESOLUTIONS:
        errors = computeErrors(M, XI_ETA_POINTS)

        pieces = []
        for label, _, _ in XI_ETA_POINTS:
            err = errors[label]
            rms = numpy.sqrt((err ** 2).mean())
            rmsByLabel[label][0].append(M)
            rmsByLabel[label][1].append(rms)
            pieces.append(err)

        pooled = numpy.concatenate(pieces)
        pooledRms = numpy.sqrt((pooled ** 2).mean())
        rmsByLabel[POOLED_LABEL][0].append(M)
        rmsByLabel[POOLED_LABEL][1].append(pooledRms)

        summary = ', '.join(f'{label.split("(")[0].strip()}={rmsByLabel[label][1][-1]:.3e}'
                             for label, _, _ in XI_ETA_POINTS)
        print(f'M={M:3d} (ncells={6 * M * M:6d}): {summary}, pooled={pooledRms:.3e}')

    alphas = plotAndFitAlpha(rmsByLabel)

    print('\nfitted error ~ M^{-alpha}:')
    for label, alpha in alphas.items():
        print(f'  {label:35s} alpha = {alpha:.3f}')

    return alphas


if __name__ == '__main__':
    main()
