#!/usr/bin/env python
"""
Maps the reconstruction-error convergence exponent alpha(xi, eta) -- as
fit in scalar_potential_reconstruction_convergence.py, but now over a fine
grid of parametric target locations covering an entire cell, not just the 4
sampled points there -- to see the actual SHAPE of the superconvergent
region for mint.VectorInterp's Whitney-1-form reconstruction of v = grad(chi).

Motivation (see the conversation this came out of): the earlier 4-point
script found alpha ~ 1.87 (2nd order) at the cell centre (xi, eta) = (0.5,
0.5), dropping to alpha ~ 0.98 (1st order) at 3 other, off-centre locations.
That's consistent with known Whitney/Nedelec edge-element theory: the
lowest-order edge element is nominally only 1st-order accurate for a
general smooth field, with 2nd-order SUPERCONVERGENCE at specific points
where the leading first-order truncation term cancels by symmetry -- the
edge-element analogue of Raviart-Thomas flux superconvergence at edge
midpoints, or Gauss-point superconvergence in Lagrange FEM. The open
question a 4-point sample can't answer: is (0.5, 0.5) the ONLY such point,
or part of a broader symmetric locus (e.g. the whole line xi=0.5, or
eta=0.5, or the diagonal)? This script fits alpha at every point of a fine
(xi, eta) grid and plots it as a heatmap to find out directly, rather than
guessing from theory.

Field, grid, and edge-data machinery are all reused UNMODIFIED from
scalar_potential_reconstruction_convergence.py (v = grad(chi), chi =
sin(lambda)*cos(theta)**2 -- exact edge data, no sympy needed, matching
mint/tests/test_vector_interp_cubedsphere.py).

Performance note: unlike the two convergence scripts, this needs MANY
target sets per resolution M (one per (xi, eta) grid point), so
mint.VectorInterp's locator is built ONCE per M and reused across many
findPoints calls -- verified directly (bit-identical results, both
numBad=0) against the fresh-locator-per-call pattern the other two scripts
use, before relying on it here for speed.

Usage: python scripts/scalar_potential_reconstruction_alpha_heatmap.py
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
from generate_cubedsphere_grid import buildCubedSphereGrid  # noqa: E402
from scalar_potential_reconstruction_convergence import (  # noqa: E402
    buildExactEdgeData, exactGradient, bilinearCellTargets, fitAlpha)

# ---------------------------------------------------------------------------
# Resolutions to fit alpha over, and the (xi, eta) grid to sample. Kept away
# from 0/1 by the same margin reasoning as the convergence scripts (see
# bilinearCellTargets' docstring): both the FindCell edge/corner tie-break
# issue, and this script's own planar-bilinear approximation of mint's
# (spherical-patch) cell geometry needing room to still land inside the
# right cell.
# ---------------------------------------------------------------------------
RESOLUTIONS = (8, 16, 32, 64)
N_GRID = 17  # xi/eta samples per axis -> N_GRID**2 locations total
MARGIN = 0.05
XI_VALS = numpy.linspace(MARGIN, 1. - MARGIN, N_GRID)
ETA_VALS = numpy.linspace(MARGIN, 1. - MARGIN, N_GRID)


def rmsErrorGrid(M):
    """
    For this M, build the grid/edge-data/locator ONCE, then for every
    (xi, eta) in the XI_VALS x ETA_VALS grid, reconstruct the vector at
    that location in every cell and return its RMS error.

    :returns: (N_GRID, N_GRID) array, rms[i, j] for (XI_VALS[i], ETA_VALS[j])
    """
    grid = buildCubedSphereGrid(M)
    data = buildExactEdgeData(grid)
    points = grid.getPoints()

    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)

    rms = numpy.zeros((N_GRID, N_GRID))
    for i, xi in enumerate(XI_VALS):
        for j, eta in enumerate(ETA_VALS):
            targetsLonLat = bilinearCellTargets(points, xi, eta)
            targets = numpy.zeros((targetsLonLat.shape[0], 3))
            targets[:, :2] = targetsLonLat

            numBad = vi.findPoints(targets, tol2=1.e-8)
            assert numBad == 0, (f'M={M}: {numBad} target points at '
                                  f'(xi={xi:.3f}, eta={eta:.3f}) fell outside their own cell')

            numeric = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
            exact = exactGradient(targetsLonLat[:, 0], targetsLonLat[:, 1])
            err = numpy.linalg.norm(numeric - exact, axis=-1)
            rms[i, j] = numpy.sqrt((err ** 2).mean())
    return rms


def main():
    # errsByM[k] = (N_GRID, N_GRID) RMS-error grid at RESOLUTIONS[k]
    errsByM = []
    for M in RESOLUTIONS:
        rms = rmsErrorGrid(M)
        errsByM.append(rms)
        print(f'M={M:3d}: rms error range over (xi, eta) grid = '
              f'[{rms.min():.3e}, {rms.max():.3e}]')
    errsByM = numpy.stack(errsByM, axis=0)  # (nM, N_GRID, N_GRID)

    # Fit alpha at every (xi, eta) grid point independently.
    alpha = numpy.zeros((N_GRID, N_GRID))
    for i in range(N_GRID):
        for j in range(N_GRID):
            a, _ = fitAlpha(RESOLUTIONS, errsByM[:, i, j])
            alpha[i, j] = a

    iC = int(numpy.argmin(numpy.abs(XI_VALS - 0.5)))
    jC = int(numpy.argmin(numpy.abs(ETA_VALS - 0.5)))
    print(f'\nalpha at (xi, eta) = ({XI_VALS[iC]:.3f}, {ETA_VALS[jC]:.3f}) [nearest to centre]: '
          f'{alpha[iC, jC]:.3f}')
    iMax, jMax = numpy.unravel_index(numpy.argmax(alpha), alpha.shape)
    print(f'max alpha = {alpha[iMax, jMax]:.3f} at (xi, eta) = '
          f'({XI_VALS[iMax]:.3f}, {ETA_VALS[jMax]:.3f})')
    print(f'median alpha over the whole grid = {numpy.median(alpha):.3f}')
    highCount = int((alpha > 1.5).sum())
    print(f'grid points with alpha > 1.5 (i.e. closer to 2nd than 1st order): '
          f'{highCount} / {N_GRID * N_GRID}')

    plotHeatmap(alpha)
    return alpha


def plotHeatmap(alpha):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # alpha[i, j] is at (XI_VALS[i], ETA_VALS[j]); pcolormesh expects
    # Z[row, col] against (x[col], y[row]) when X/Y are the mesh edges, so
    # transpose to put xi on the x-axis, eta on the y-axis as usual.
    mesh = ax.pcolormesh(XI_VALS, ETA_VALS, alpha.T, shading='nearest',
                          cmap='viridis', vmin=0.5, vmax=2.0)
    fig.colorbar(mesh, ax=ax, label='fitted alpha (error ~ M^-alpha)')

    cs = ax.contour(XI_VALS, ETA_VALS, alpha.T, levels=[1.0, 1.5],
                     colors=['white', 'red'], linewidths=1.2)
    ax.clabel(cs, inline=True, fontsize=8, fmt='%.1f')

    ax.plot([0.5], [0.5], marker='*', markersize=18, color='red',
             markeredgecolor='black', label='(xi, eta) = (0.5, 0.5)')
    ax.set_xlabel('xi')
    ax.set_ylabel('eta')
    ax.set_aspect('equal')
    ax.set_title('Convergence exponent alpha(xi, eta)\n'
                  r'$\mathbf{v} = \nabla\chi$, $\chi=\sin\lambda\,\cos^2\theta$, '
                  f'M = {RESOLUTIONS}')
    ax.legend(loc='upper right', fontsize='small')
    fig.tight_layout()

    outfile = HERE / 'scalar_potential_reconstruction_alpha_heatmap.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
