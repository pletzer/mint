#!/usr/bin/env python
"""
Convergence of the 2nd-order (k=2) mimetic/Whitney edge element -- see
scripts/whitney2.py -- on a plain Cartesian x,y grid, as a fast first check
before trying the same thing on the cubed-sphere
(whitney2_convergence_cubedsphere.py). No mint dependency at all: this
script builds its own uniform grid directly.

Test field: u = mimetic/Whitney 2 self-test:

    u = curl(psi * zhat) = (dpsi/dy, -dpsi/dx),   psi(x, y) = sin(KX*x)*sin(KY*y)

(a genuinely rotational, non-gradient field -- circulations have to
actually be integrated, not just differenced from a potential -- matching
the convention used throughout mint/scripts/pole_asymmetric_convergence.py
and vector_potential_reconstruction_convergence.py).

Per cell, the reconstructed line integral along a fixed, deliberately
off-centre chord in (xi, eta) (whitney2.DEFAULT_CHORD -- NOT the
corner-to-corner diagonal, which passes through the cell centroid, an
isolated superconvergent point that masks the difference between k=1 and
k=2, see whitney2.py's module docstring) is compared against the true
line integral along that same chord (whitney2.exactChordLineIntegral),
for all three slope-DOF strategies ('order1' i.e. plain k=1 Whitney,
'exact', 'split' -- see whitney2.buildCellDOFs), and the RMS error over
all cells is fit to error ~ N^-alpha as the grid is refined.

Expectation, confirmed both symbolically (Taylor-expanding a generic
smooth field around a cell corner) and empirically here: 'order1' has a
genuine leading O(h^2) error term; for 'exact'/'split' that h^2 term
cancels EXACTLY (not just approximately -- confirmed symbolically, for
any smooth field), leaving a leading O(h^3) error -- i.e. the slope DOF
buys exactly ONE extra order, not two (an earlier, sloppier expectation
of O(h^2) -> O(h^4) in this codebase's history was wrong).

That one-order gain is real but easy to miss in a narrow-range fit: for
THIS field, order1's h^2 coefficient happens to be small relative to its
own h^3 coefficient (which picks up extra factors from second
derivatives/wavenumbers-squared), so at moderate N its error is still
riding on that h^3 term -- 'order1' and 'exact'/'split' look deceptively
similar (fitted alpha ~2.9-3.0 for ALL three modes) until N is pushed high
enough (empirically, past ~64 for this field) that order1's h^2 term
finally takes over and its LOCAL (pairwise, log2(err[N]/err[2N])) order
visibly bends down from ~2.9 towards 2, while 'exact'/'split' lock onto
exactly 3.00 and stay there. RESOLUTIONS below is deliberately extended
into that regime to make the bend visible rather than just fitting a
single alpha over a range that is still pre-asymptotic for order1.

Usage: python scripts/whitney2_convergence_cartesian.py
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

LX, LY = 1.0, 1.0
KX, KY = 2.0 * numpy.pi, numpy.pi


def vectorField(pts):
    """u = dpsi/dy, w = -dpsi/dx for psi(x, y) = sin(KX*x)*sin(KY*y)."""
    x, y = pts[..., 0], pts[..., 1]
    u = KY * numpy.sin(KX * x) * numpy.cos(KY * y)
    w = -KX * numpy.cos(KX * x) * numpy.sin(KY * y)
    return numpy.stack([u, w], axis=-1)


def linearVectorField(pts):
    """
    psi = x*y => u = dpsi/dy = x, w = -dpsi/dx = -y: affine in (x, y), so
    along ANY straight edge (any orientation) the tangential component is
    exactly linear in arclength fraction s -- used by selfTestLinearField
    to check the k=2 element reproduces such a field's diagonal line
    integral to machine precision.
    """
    x, y = pts[..., 0], pts[..., 1]
    return numpy.stack([x, -y], axis=-1)


def buildGrid(N):
    """
    Uniform N x N grid on [0, LX] x [0, LY]; returns (N*N, 4, 2) corners
    in the 0,1,2,3 convention (see whitney2.buildCellDOFs docstring).
    """
    xs = numpy.linspace(0.0, LX, N + 1)
    ys = numpy.linspace(0.0, LY, N + 1)
    corners = numpy.empty((N * N, 4, 2))
    icell = 0
    for i in range(N):
        for j in range(N):
            corners[icell, 0] = (xs[i], ys[j])
            corners[icell, 1] = (xs[i + 1], ys[j])
            corners[icell, 2] = (xs[i + 1], ys[j + 1])
            corners[icell, 3] = (xs[i], ys[j + 1])
            icell += 1
    return corners


def selfTestLinearField():
    """
    Regression check for the sign subtlety flagged in whitney2.py's module
    docstring: for a field that is affine in (x, y) (hence exactly linear
    along every straight edge), both 'exact' and 'split' should reproduce
    the diagonal line integral to near machine precision, while 'order1'
    (no slope term) should NOT. Uses a genuine (non-rectangular,
    non-axis-aligned) parallelogram cell: for an affine bilinear map the
    k=2 tensor-product ansatz is exactly rich enough to represent an
    affine field with zero residual (verified analytically -- a general,
    non-parallelogram quad would leave a small residual from the map's own
    bilinear curvature, which is not what this check is targeting).
    """
    corners = numpy.array([[0.3, 0.7], [1.1, 0.9], [1.3, 1.7], [0.5, 1.5]])
    exact = exactChordLineIntegral(corners, linearVectorField)
    for mode, tol in (('exact', 1e-9), ('split', 1e-9), ('order1', None)):
        a, b = buildCellDOFs(corners, linearVectorField, mode=mode)
        recon = chordLineIntegral(a, b)
        err = abs(recon - exact)
        print(f'  self-test [{mode:6s}]: reconstructed={recon:.10f}  exact={exact:.10f}  err={err:.2e}')
        if tol is not None:
            assert err < tol, f'self-test FAILED for mode={mode}: err={err:.2e} >= {tol:.2e}'
    print('  self-test passed (exact/split reproduce an affine field to machine precision;'
          ' order1 does not, as expected)\n')


RESOLUTIONS = (4, 8, 16, 32, 64, 128, 256)
MODES = ('order1', 'exact', 'split')


def computeRmsErrors(N):
    corners = buildGrid(N)
    ncells = corners.shape[0]
    errs = {mode: numpy.empty(ncells) for mode in MODES}
    for icell in range(ncells):
        c = corners[icell]
        exact = exactChordLineIntegral(c, vectorField)
        for mode in MODES:
            a, b = buildCellDOFs(c, vectorField, mode=mode)
            recon = chordLineIntegral(a, b)
            errs[mode][icell] = recon - exact
    return {mode: float(numpy.sqrt(numpy.mean(errs[mode] ** 2))) for mode in MODES}


def fitAlpha(Ns, errs):
    """Least-squares fit of log(err) = -alpha*log(N) + c; returns (alpha, prefactor)."""
    Ns = numpy.asarray(Ns, dtype=numpy.float64)
    errs = numpy.asarray(errs, dtype=numpy.float64)
    slope, intercept = numpy.polyfit(numpy.log(Ns), numpy.log(errs), 1)
    return -slope, numpy.exp(intercept)


def localOrders(rms_mode):
    """Pairwise (local) order log2(err[N]/err[2N]) between successive
    RESOLUTIONS -- more diagnostic than a single global fit when the
    global fit straddles a pre-asymptotic/asymptotic crossover (see module
    docstring: this is exactly what happens for 'order1' here)."""
    return [numpy.log2(rms_mode[i] / rms_mode[i + 1]) for i in range(len(rms_mode) - 1)]


def main():
    print('self-test on a single parallelogram cell with an affine field:')
    selfTestLinearField()

    rms = {mode: [] for mode in MODES}
    for N in RESOLUTIONS:
        r = computeRmsErrors(N)
        for mode in MODES:
            rms[mode].append(r[mode])
        print(f'N={N:4d} (ncells={N * N:6d}): ' + ', '.join(f'{m}={r[m]:.3e}' for m in MODES))

    print('\nlocal (pairwise) order = log2(err[N]/err[2N]) -- watch order1 bend '
          'from ~3 towards 2 as N grows, while exact/split lock onto 3:')
    for mode in MODES:
        orders = ', '.join(f'{o:.2f}' for o in localOrders(rms[mode]))
        print(f'  {mode:8s}: {orders}')

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {'order1': 'tab:blue', 'exact': 'tab:green', 'split': 'tab:red'}
    Ns = numpy.asarray(RESOLUTIONS, dtype=numpy.float64)
    print()
    for mode in MODES:
        alpha, prefactor = fitAlpha(RESOLUTIONS, rms[mode])
        ax.loglog(Ns, rms[mode], 'o-', color=colors[mode], label=f'{mode} (alpha={alpha:.2f})')
        ax.loglog(Ns, prefactor * Ns ** (-alpha), '--', color=colors[mode], alpha=0.5, linewidth=1)
        print(f'{mode:8s}: fitted error ~ N^-{alpha:.3f}')

    ax.set_xlabel('N  (cells per side; N^2 cells total)')
    ax.set_ylabel('RMS chord line-integral error')
    ax.set_title('2nd-order Whitney convergence, Cartesian grid\n'
                  r'$u=\partial_y\psi,\ w=-\partial_x\psi,\ \psi=\sin(2\pi x)\sin(\pi y)$'
                  '\n(global fit over the whole range; see printed LOCAL orders for the'
                  '\norder1 h^2-vs-h^3 crossover this fit alone hides)')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'whitney2_convergence_cartesian.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
