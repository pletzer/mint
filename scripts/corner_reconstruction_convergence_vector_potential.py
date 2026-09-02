#!/usr/bin/env python
"""
Reconstructs u = curl(psi * r_hat) AT cubed-sphere GRID CORNER (vertex)
points, via a "+" stencil: two short mint.PolylineIntegral line integrals
centred exactly at the corner, one running west-east and one south-north,
in PURE (lon, lat) coordinates (not along any particular cell's own edges),
each spanning one nominal cell-width. Dividing each branch's line integral
by its own COORDINATE length (plain degrees, matching mint's raw-degree
convention -- see below) recovers an estimate of the field's raw
(d/dlon_deg, d/dlat_deg) component there.

Why this should be O(h^2) despite pointwise Whitney reconstruction only
being O(h) away from a cell's own centroid (see
scalar_potential_reconstruction_alpha_heatmap.py): worked out in the
conversation this came out of, short version below.

By the "exact midpoint" identity (see vector_potential_reconstruction_
convergence.py's sibling scripts' docstrings), each HALF of the west-east
branch -- a straight path from the corner to a point delta away, entirely
within ONE cell -- has its line integral EXACTLY equal to v_h evaluated at
that half-branch's own parametric midpoint, times delta. That midpoint
sits at parametric offset ~+delta/2's worth of distance from the WEST
cell's own centroid, and ~-delta/2's worth from the EAST cell's centroid
(the two cells are mirror images of each other about the shared corner).
Writing out the Taylor error of the Whitney reconstruction at each of those
two points (see scalar_potential_reconstruction_alpha_heatmap.py's finding
that this error's leading term is LINEAR in the offset from the cell's own
centroid):

  [v_h(+off) - v_true(+off)]  +  [v_h(-off) - v_true(-off)]  =  O(h^2)

even though each individual bracket is only O(h) -- the leading terms are
equal in magnitude and opposite in sign (same underlying smooth field,
offsets of opposite sign), so they cancel when the two half-branches are
combined into one full branch's average. This is the corner-stencil analog
of the diagonal-path result: instead of landing on the ONE superconvergent
point inside a single cell, it cancels the leading error algebraically
between two neighbouring cells.

Units note (a repeat of the vector_potential_reconstruction_convergence.py
DEG2RAD gotcha, in a new guise): mint's line integral is computed using the
grid's raw, DEGREE-valued point coordinates as-is -- no spherical metric
correction. So "dividing by branch length" must use the branch's plain
COORDINATE length (e.g. 2*halfWidthDeg for the west-east branch, NOT a
metric-corrected physical arc length with a cos(lat) factor), to stay
consistent with exactVector's own raw-degree-basis convention (see
vector_potential_reconstruction_convergence.py's exactVector docstring).

Test field: u = curl(psi * r_hat), psi = cos(theta)(1+sin(theta))cos(lambda)
-- the same "vector potential" field used throughout this codebase. Edge
data is EXACT (sympy, no quadrature) via vector_potential_reconstruction_
convergence.buildExactEdgeData -- same cost/RESOLUTIONS constraint as that
script (~60ms/edge measured there; capped at (4, 8, 16, 32); run this in
the background, expect ~10-35 minutes).

Four test corners per M, chosen to probe increasingly singular grid
topology, located directly from grid.getPoints() (cell (i, j)'s own
corner-c is exactly logical vertex (i, j) shifted by c -- see
cubedSphereGridPoints/cellIndex below). All four sit at FIXED (lon, lat)
independent of M (exact cube-geometric points, not just M-converging
approximations), verified directly before running the expensive sympy
edge data:

  - 'interior'        (0, 0):        (lon,lat) = (0, 0) on '+X' -- shared
                        by 4 cells, all on the SAME panel.
  - '2-panel edge'     (-45, 0):      the seam between '+X' and '-Y' --
                        its west branch crosses into the neighbouring panel.
  - '3-panel junction' (45, 35.264):  the cube corner where '+X', '+Y' and
                        '+Z' all meet (the (1,1,1)/sqrt(3) direction) --
                        only 3 cells share this vertex, not 4 or 2.
  - 'near dateline'    (180, 0):      centre of the '-X' panel, sitting
                        EXACTLY on the +-180 branch cut -- the west-east
                        branch's east endpoint is at 180+delta > 180,
                        past the wraparound boundary. Verified separately
                        (not asserted blindly) that mint.PolylineIntegral
                        with periodX=360 handles this identically whether
                        fed as literal 180+delta or pre-wrapped to
                        -180+delta (bit-identical flux either way) -- so
                        no special-casing needed, the raw unwrapped value
                        is used directly below.

Usage: python scripts/corner_reconstruction_convergence_vector_potential.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import PolylineIntegral, CELL_BY_CELL_DATA

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from generate_cubedsphere_grid import buildCubedSphereGrid  # noqa: E402
from vector_potential_reconstruction_convergence import buildExactEdgeData, exactVector  # noqa: E402
from scalar_potential_reconstruction_convergence import fitAlpha  # noqa: E402 (generic, field-agnostic helper)

RESOLUTIONS = (4, 8, 16, 32)  # not higher: exactLineIntegral is per-edge sympy, see module docstring

# (label, (faceIdx, i, j, cornerIdx) as a function of M) -- faceIdx indexes
# generate_cubedsphere_grid._faceFrames ('+X'=0, '-X'=1, '+Y'=2, '-Y'=3,
# '+Z'=4, '-Z'=5); cell (i,j)'s corner-cornerIdx is vertex (i,j) shifted by
# cornerIdx per cubedSphereGridPoints' corner order (0=(i,j), 1=(i+1,j),
# 2=(i+1,j+1), 3=(i,j+1)).
TEST_CORNERS = [
    ('interior',        lambda M: (0, M // 2, M // 2, 0)),
    ('2-panel edge',     lambda M: (0, 0, M // 2, 0)),
    ('3-panel junction', lambda M: (0, M - 1, M - 1, 2)),
    ('near dateline',    lambda M: (1, M // 2, M // 2, 0)),
]


def cellIndex(faceIdx, i, j, M):
    """Cell (i, j)'s flat index on face faceIdx, matching cubedSphereGridPoints's build order."""
    return faceIdx * M * M + i * M + j


def branchIntegral(grid, data, p0, p1):
    """mint.PolylineIntegral's flux of `data` along the straight path p0->p1 (lon,lat degrees)."""
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
    xyz = numpy.array([(p0[0], p0[1], 0.), (p1[0], p1[1], 0.)])
    pli.computeWeights(xyz, counterclock=False)
    return pli.getIntegral(data, placement=CELL_BY_CELL_DATA)


def cornerEstimate(grid, data, lon0, lat0, halfWidthDeg):
    """
    The "+" stencil: (u_est, v_est) at (lon0, lat0), in mint's raw
    (d/dlon_deg, d/dlat_deg) coordinate-basis convention (see module
    docstring's units note).
    """
    we = branchIntegral(grid, data, (lon0 - halfWidthDeg, lat0), (lon0 + halfWidthDeg, lat0))
    sn = branchIntegral(grid, data, (lon0, lat0 - halfWidthDeg), (lon0, lat0 + halfWidthDeg))
    return we / (2. * halfWidthDeg), sn / (2. * halfWidthDeg)


def main():
    errsByLabel = {label: ([], []) for label, _ in TEST_CORNERS}

    for M in RESOLUTIONS:
        print(f'M={M}: building exact (sympy) edge data for {6 * M * M} cells '
              f'({6 * M * M * 4} edges, slow)...', flush=True)
        grid = buildCubedSphereGrid(M)
        data = buildExactEdgeData(grid)
        points = grid.getPoints()

        halfWidthDeg = 0.5 * (90.0 / M)  # one nominal cell-width, centred at the corner

        for label, cornerFn in TEST_CORNERS:
            faceIdx, I, J, cornerIdx = cornerFn(M)
            lon0, lat0 = points[cellIndex(faceIdx, I, J, M), cornerIdx, :2]

            u_est, v_est = cornerEstimate(grid, data, lon0, lat0, halfWidthDeg)
            exact = exactVector(lon0, lat0)
            err = numpy.hypot(u_est - exact[0], v_est - exact[1])

            errsByLabel[label][0].append(M)
            errsByLabel[label][1].append(err)
            print(f'  M={M:3d} {label:18s} corner=({lon0:8.3f},{lat0:7.3f})  '
                  f'u_est={u_est: .5f} v_est={v_est: .5f}  '
                  f'exact=({exact[0]: .5f},{exact[1]: .5f})  err={err:.3e}', flush=True)

    print()
    alphas = {}
    for label, (Ms, errs) in errsByLabel.items():
        alpha, prefactor = fitAlpha(Ms, errs)
        alphas[label] = alpha
        print(f'{label:18s}: fitted error ~ M^-alpha, alpha = {alpha:.3f}')

    plotResult(errsByLabel, alphas)
    return alphas


def plotResult(errsByLabel, alphas):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    colors = plt.cm.tab10.colors

    for i, (label, (Ms, errs)) in enumerate(errsByLabel.items()):
        Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
        errs_arr = numpy.asarray(errs, dtype=numpy.float64)
        alpha, prefactor = fitAlpha(Ms, errs)
        color = colors[i % len(colors)]
        ax.loglog(Ms_arr, errs_arr, 'o-', color=color, linewidth=2,
                   label=f'{label}  (alpha={alpha:.2f})')
        ax.loglog(Ms_arr, prefactor * Ms_arr ** (-alpha), '--', color=color, linewidth=1, alpha=0.6)

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel(r'$|\mathbf{v}_{stencil} - \mathbf{v}_{exact}|$ at corner')
    ax.set_title('"+" stencil corner reconstruction error vs resolution\n'
                  r'$\mathbf{u} = \nabla\times(\psi\,\hat{r})$, '
                  r'$\psi=\cos\theta(1+\sin\theta)\cos\lambda$')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'corner_reconstruction_convergence_vector_potential.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
