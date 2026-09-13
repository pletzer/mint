#!/usr/bin/env python
"""
Convergence of the 2nd-order (k=2) mimetic/Whitney edge element (see
scripts/whitney2.py and whitney2_convergence_cartesian.py's Cartesian
version, run FIRST) on cubed-sphere grids, using the same pole-asymmetric
vector-potential field u = curl(psi*r_hat),
psi = cos(theta)*(1+sin(theta))*cos(lambda), already used throughout this
codebase (pole_asymmetric_convergence.py,
vector_potential_reconstruction_convergence.py).

Unlike vector_potential_reconstruction_convergence.py (which builds edge
data via a fresh sympy symbolic integration per edge, exactLineIntegral,
~60ms/edge for this field), this script gets circulations and slopes via
plain Gauss-Legendre quadrature directly on the closed-form vectorField
(reused unmodified from vector_potential_reconstruction_convergence.exactVector,
which already carries the DEG2RAD raw-coordinate-basis correction that
module's docstring derives) -- much cheaper (no sympy at all), so this can
go to higher M within a reasonable runtime.

Rather than a single cross-panel path (as in
line_integral_convergence_vector_potential.py, which needs mint's own
locator to find which cell each point of an arbitrary path falls in),
this measures the RMS line-integral error over EVERY cell's own fixed,
off-centre (xi, eta) chord (whitney2.DEFAULT_CHORD -- NOT the corner
diagonal, see whitney2.py's module docstring for why) -- avoiding any
dependency on mint's locator/point-location API, and directly comparable
in style to vector_potential_reconstruction_convergence.py's per-cell RMS
approach.

Expectation, per whitney2_convergence_cartesian.py's own result: don't
expect a dramatic O(h^2) -> O(h^4) jump even where this DOES improve --
the honest theoretical gain from adding the slope DOF is one order
(O(h^2) -> O(h^3), see that script's module docstring), and even that can
be hidden by a pre-asymptotic crossover at the resolutions tested.

PANEL RESTRICTION: buildCubedSphereGrid lays out cells face-by-face, as 6
CONTIGUOUS blocks of M*M cells each, in the fixed order
['+X', '-X', '+Y', '-Y', '+Z', '-Z'] (see cubedSphereGridPoints's loop
over _faceFrames -- same order as generate_cubedsphere_grid._CUBE_FACE_NAMES,
duplicated here rather than importing that private name). '+Z'/'-Z' are
the two panels containing a pole; '+X','-X','+Y','-Y' are purely
equatorial (their own corners never exceed ~35.26 deg latitude, the
cube-corner projection angle) and never see the extreme lon/lat cell
distortion a pole forces on its own panel. selectPanelCells below slices
out one such block directly (cheap, exact) and cross-checks the slice
against generate_cubedsphere_grid.cubeFaceOf on every cell's own corner-0
(catches an ordering-assumption bug immediately rather than silently
testing the wrong cells).

This lets computeRmsErrors optionally restrict to '+X' alone, to test the
hypothesis that the missing improvement on the FULL sphere (previous run
of this script) is a lon/lat coordinate-distortion effect concentrated at
the poles, rather than the geometric (non-parallelogram, no interior
bubble DOFs) limitation flagged when the k=2 idea was first proposed.

RESULT (confirms the pole hypothesis): on the whole sphere, order1 AND
exact/split all converge to the SAME ~O(h^2) ceiling (local order -> 1.9-2.0
for every mode, even at M=64) -- the pole panels cap k=2's benefit just as
much as k=1's, consistent with the missing-bubble-DOF limitation being
concentrated right at the singular corners. Restricted to the '+X' panel
alone (no pole, 6x cheaper so pushed to M=256), the SAME crossover already
found on the Cartesian grid reappears: order1's local order visibly bends
down from ~3 towards 2 (2.96, 2.79, 2.43, 2.13 as M: 16->32->64->128->256),
while exact/split lock onto exactly 3.00 from M=32 onward. So the k=2
element's genuine one-order gain IS present on well-behaved (non-polar)
cubed-sphere cells -- it was invisible on the full-sphere test purely
because the pole panels (2 of 6) impose their own, unrelated O(h^2) ceiling
that the RMS-over-the-whole-sphere norm is dominated by as h -> 0 (a
global error norm is limited by its worst-behaved subregion).

Usage: python scripts/whitney2_convergence_cubedsphere.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from whitney2 import buildCellDOFs, chordLineIntegral, exactChordLineIntegral  # noqa: E402
from generate_cubedsphere_grid import buildCubedSphereGrid, cubeFaceOf  # noqa: E402
from vector_potential_reconstruction_convergence import exactVector  # noqa: E402

# same order cubedSphereGridPoints builds panels in (see module docstring)
_PANEL_ORDER = ('+X', '-X', '+Y', '-Y', '+Z', '-Z')


def vectorField(pts):
    """(..., 2) (lon_deg, lat_deg) -> (..., 2) raw-coordinate-basis vector
    (u_lam, u_theta), reusing the already-verified exactVector."""
    v3 = exactVector(pts[..., 0], pts[..., 1])
    return v3[..., :2]


MODES = ('order1', 'exact', 'split')


def selectPanelCells(points, M, panel):
    """
    points: (6*M*M, 4, 2) as returned by grid.getPoints()[..., :2].
    panel: one of _PANEL_ORDER, or None for all cells (no-op).

    Slices out the ncells==M*M contiguous block belonging to `panel`
    (cheap/exact -- see module docstring for the ordering this relies on),
    then cross-checks a sample of cells against cubeFaceOf so an ordering
    bug would fail loudly instead of silently testing the wrong panel.
    """
    if panel is None:
        return points
    ncells_per_panel = M * M
    ipanel = _PANEL_ORDER.index(panel)
    block = points[ipanel * ncells_per_panel:(ipanel + 1) * ncells_per_panel]

    sampleIdx = numpy.linspace(0, block.shape[0] - 1, min(8, block.shape[0])).astype(int)
    for i in sampleIdx:
        # cell CENTROID, not a corner -- a cell's own corner-0 grid point can
        # be a genuine cube VERTEX (equidistant from 3 faces, e.g. the '+X'
        # panel's own (i=0,j=0) corner ties '+X'/'-Y'/'-Z'), where cubeFaceOf's
        # argmax is a coin flip on floating-point noise; the centroid is
        # always safely interior to the panel.
        centroidLon, centroidLat = block[i, :, 0].mean(), block[i, :, 1].mean()
        face = cubeFaceOf(centroidLon, centroidLat)
        assert face == panel, (
            f'panel slice mismatch: cell {i} of the assumed {panel!r} block is '
            f'actually on face {face!r} -- the panel ordering assumption in this '
            f'module\'s docstring no longer holds, fix selectPanelCells')
    return block


def computeRmsErrors(M, panel=None):
    grid = buildCubedSphereGrid(M)
    points = selectPanelCells(grid.getPoints()[..., :2], M, panel)  # (ncells, 4, 2): lon, lat in degrees
    ncells = points.shape[0]
    errs = {mode: numpy.empty(ncells) for mode in MODES}
    for icell in range(ncells):
        c = points[icell]
        exact = exactChordLineIntegral(c, vectorField)
        for mode in MODES:
            a, b = buildCellDOFs(c, vectorField, mode=mode)
            recon = chordLineIntegral(a, b)
            errs[mode][icell] = recon - exact
    return {mode: float(numpy.sqrt(numpy.mean(errs[mode] ** 2))) for mode in MODES}


def fitAlpha(Ms, errs):
    """Least-squares fit of log(err) = -alpha*log(M) + c; returns (alpha, prefactor)."""
    Ms = numpy.asarray(Ms, dtype=numpy.float64)
    errs = numpy.asarray(errs, dtype=numpy.float64)
    slope, intercept = numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)
    return -slope, numpy.exp(intercept)


def localOrders(rms_mode):
    """Pairwise (local) order log2(err[M]/err[2M]) between successive
    RESOLUTIONS -- see whitney2_convergence_cartesian.py's localOrders for
    why this is more diagnostic than a single global fit."""
    return [numpy.log2(rms_mode[i] / rms_mode[i + 1]) for i in range(len(rms_mode) - 1)]


SCOPES = (
    # (label, panel, ncells-per-M factor, resolutions)
    # Whole-sphere kept to M<=64 (6*M^2 cells, real cost); the equatorial
    # panel is M^2 cells (6x cheaper), pushed much further (M=256) since
    # that's what it takes to see order1 visibly bend from ~3 towards 2
    # while exact/split hold at 3 -- the SAME crossover already found on
    # the Cartesian grid, now isolated to the pole-free part of the sphere.
    ('all panels (whole sphere, incl. poles)', None, 6, (4, 8, 16, 32, 64)),
    ("'+X' panel only (equatorial, no pole)", '+X', 1, (4, 8, 16, 32, 64, 128, 256)),
)


def main():
    results = {}
    for label, panel, ncellsFactor, resolutions in SCOPES:
        print(f'--- {label} ---')
        rms = {mode: [] for mode in MODES}
        for M in resolutions:
            r = computeRmsErrors(M, panel=panel)
            for mode in MODES:
                rms[mode].append(r[mode])
            print(f'M={M:3d} (ncells={ncellsFactor * M * M:6d}): '
                  + ', '.join(f'{m}={r[m]:.3e}' for m in MODES))
        print('local (pairwise) order = log2(err[M]/err[2M]):')
        for mode in MODES:
            orders = ', '.join(f'{o:.2f}' for o in localOrders(rms[mode]))
            print(f'  {mode:8s}: {orders}')
        print()
        results[label] = (resolutions, rms)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    colors = {'order1': 'tab:blue', 'exact': 'tab:green', 'split': 'tab:red'}
    linestyles = {SCOPES[0][0]: '--', SCOPES[1][0]: '-'}
    markers = {SCOPES[0][0]: 's', SCOPES[1][0]: 'o'}
    for label, (resolutions, rms) in results.items():
        Ms = numpy.asarray(resolutions, dtype=numpy.float64)
        for mode in MODES:
            alpha, _ = fitAlpha(resolutions, rms[mode])
            ax.loglog(Ms, rms[mode], linestyles[label], marker=markers[label],
                      color=colors[mode], alpha=0.85 if label == SCOPES[1][0] else 0.5,
                      label=f'{mode}, {label} (alpha={alpha:.2f})')
            print(f'{mode:8s} [{label}]: fitted error ~ M^-{alpha:.3f}')

    ax.set_xlabel('M  (cells per cubed-sphere panel edge)')
    ax.set_ylabel('RMS chord line-integral error')
    ax.set_title('2nd-order Whitney convergence, cubed-sphere\n'
                  'whole sphere (dashed) vs single equatorial panel (solid)\n'
                  r'$\mathbf{u}=\nabla\times(\psi\,\hat r)$, '
                  r'$\psi=\cos\theta(1+\sin\theta)\cos\lambda$')
    ax.legend(fontsize='x-small', loc='best')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'whitney2_convergence_cubedsphere.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
