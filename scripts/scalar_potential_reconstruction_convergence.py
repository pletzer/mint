#!/usr/bin/env python
"""
Convergence of mint.VectorInterp's vector reconstruction (from edge-integral,
1-form data) on cubed-sphere grids of increasing resolution M, for target
points placed at arbitrary parametric locations (xi, eta) inside each cell.

This is the SCALAR-potential counterpart to
scripts/vector_potential_reconstruction_convergence.py -- same grids, same
(xi, eta) target locations, same fitting/plotting machinery (reused directly
from that script, see the imports below), but a different vector field:

    v = grad(chi),   chi(lambda, theta) = sin(lambda) * cos(theta)**2

the same 0-form used throughout mint/tests/test_vector_interp_cubedsphere.py
(chi's own gradient, and the chain-rule DEG2RAD factor in exactGradient
below, are copied from that test unmodified -- see its docstring for the
derivation).

Why this comparison matters (see the vector-potential script's docstring for
the full discussion): a GRADIENT field's edge integral is exact by the
fundamental theorem of calculus -- e_i = chi(p1) - chi(p0), a plain endpoint
difference, independent of the actual shape of the path/edge between them.
So unlike the vector-potential script's u = curl(psi * r_hat) (whose edge
data comes from integrating along a flat-(lon,lat)-chord approximation of
the true, closer-to-great-circle mesh edge -- an extra O(1/M) error source
baked in before reconstruction even starts), THIS script's edge data carries
NO geometric-approximation error at all, at any M. Whatever convergence rate
comes out here is therefore the Whitney-1-form reconstruction's OWN error,
isolated -- not lumped together with a second, independent error source.
mint/tests/test_vector_interp_cubedsphere.py already found, at cell centres
only, RMS ~ 1/M^2 (alpha ~ 1.86, confirmed by direct fit) and MAX ~ 1/M
(alpha ~ 1.03, pole-adjacent cells) -- this script re-derives the RMS number
at cell centres as a sanity check (should closely reproduce ~1.86) and
extends the same measurement to the vector-potential script's other 3
off-centre (xi, eta) locations, so the two fields' alphas are directly
comparable side by side, location by location.

No sympy needed here (unlike the vector-potential script) -- endpoint
differences are cheap -- so RESOLUTIONS goes back up to include 64.

Usage: python scripts/scalar_potential_reconstruction_convergence.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import VectorInterp, CELL_BY_CELL_DATA

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from generate_cubedsphere_grid import buildCubedSphereGrid, NUM_VERTS_PER_QUAD  # noqa: E402
from vector_potential_reconstruction_convergence import bilinearCellTargets  # noqa: E402

DEG2RAD = numpy.pi / 180.0


# ---------------------------------------------------------------------------
# The scalar potential (0-form) and the gradient field it derives -- copied
# from mint/tests/test_vector_interp_cubedsphere.py (potential/exactGradient)
# unmodified, see that file's docstring for the derivation, in particular
# why exactGradient needs the extra DEG2RAD chain-rule factor (mint's
# reconstructed vector is expressed w.r.t. the grid's raw DEGREE-valued
# coordinate basis, not radians).
# ---------------------------------------------------------------------------
def potential(lonDeg, latDeg):
    """The 0-form: chi(lon, lat) = sin(lon) * cos(lat)**2."""
    lam = lonDeg * DEG2RAD
    theta = latDeg * DEG2RAD
    return numpy.sin(lam) * numpy.cos(theta) ** 2


def exactGradient(lonDeg, latDeg):
    """
    Exact grad(chi), expressed in the grid's own (lon_deg, lat_deg)
    coordinate basis (matching mint.VectorInterp's raw output). Accepts
    scalar or array lonDeg/latDeg; returns an array of shape (..., 3) with a
    trailing zero z-component.
    """
    lam = lonDeg * DEG2RAD
    theta = latDeg * DEG2RAD
    dlon = DEG2RAD * numpy.cos(lam) * numpy.cos(theta) ** 2
    dlat = -DEG2RAD * numpy.sin(lam) * numpy.sin(2. * theta)
    return numpy.stack([dlon, dlat, numpy.zeros_like(dlon)], axis=-1)


def buildExactEdgeData(grid):
    """
    Cell-by-cell edge-integrated data for dchi, EXACT for any potential and
    any edge shape by the fundamental theorem of calculus -- data[:, i0] is
    exactly chi(edge's end node) - chi(edge's start node) -- see module
    docstring. Same corner/edge ordering and sign convention as
    vector_potential_reconstruction_convergence.buildExactEdgeData:

        3-->--2
        |     |
        ^     ^
        |     |
        0-->--1
    """
    points = grid.getPoints()
    ncells = points.shape[0]
    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * (potential(points[:, i1, 0], points[:, i1, 1])
                               - potential(points[:, i0, 0], points[:, i0, 1]))
    return data


def computeErrors(M, xietaList):
    """
    Build the cubed-sphere grid and its exact edge data once for this M,
    then for each (label, xi, eta) in xietaList, locate the corresponding
    target point inside every cell and reconstruct the vector there.

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
        exact = exactGradient(targetsLonLat[:, 0], targetsLonLat[:, 1])
        errors[label] = numpy.linalg.norm(numeric - exact, axis=-1)
    return errors


# ---------------------------------------------------------------------------
# Resolutions and the set of (xi, eta) target locations -- identical to
# vector_potential_reconstruction_convergence.py, for a direct comparison.
# ---------------------------------------------------------------------------
RESOLUTIONS = (4, 8, 16, 32, 64)

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
                  r'$\mathbf{v} = \nabla\chi$, $\chi=\sin\lambda\,\cos^2\theta$')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'scalar_potential_reconstruction_convergence.png'
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
