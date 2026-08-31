#!/usr/bin/env python
"""
Convergence of mint.PolylineIntegral's line integral of v = grad(chi) along
a path that crosses multiple cubed-sphere panels, as the grid is refined.

chi(lambda, theta) = sin(lambda) * cos(theta)**2 -- the same 0-form used
throughout this codebase's cubed-sphere tests (see
mint/tests/test_vector_interp_cubedsphere.py and
scalar_potential_reconstruction_convergence.py, whose buildExactEdgeData is
reused here unmodified). Edge data is EXACT: e_i = chi(edge end) -
chi(edge start), by the fundamental theorem of calculus, no quadrature or
sympy needed at all.

The REFERENCE value is exact too, and even cheaper: for a gradient field,
the line integral along ANY path -- regardless of its shape, or how many
grid cells/panels it crosses -- is exactly chi(pathEnd) - chi(pathStart),
by the same fundamental theorem of calculus. Just two function evaluations,
independent of M.

Path: (-30, -20) -> (45, 50), reused from
scripts/pole_asymmetric_convergence.py's 'away from poles' test line,
VERIFIED below (generate_cubedsphere_grid.facesAlongPath) to cross from the
+X cube face into the +Z (north pole) face -- i.e. it genuinely crosses a
cubed-sphere panel seam, not just many cells within one panel.

Prediction (see the conversation this came out of -- follows on from
scalar_potential_reconstruction_alpha_heatmap.py, which found the pointwise
RECONSTRUCTED VECTOR is only 1st-order accurate almost everywhere in a
cell, with 2nd order isolated to the single point (xi, eta) = (0.5, 0.5)):
a LINE INTEGRAL is fundamentally a difference of the underlying
bilinearly-interpolated POTENTIAL at the path's two endpoints (mint's
lowest-order Whitney edge elements are the exact "d" of Whitney nodal (Q1)
elements), and standard nodal bilinear interpolation of a smooth function
is UNIFORMLY 2nd-order accurate -- not just at one special point. So,
unlike the pointwise-vector study, this is expected to converge at
alpha ~ 2 regardless of exactly how the path sits relative to any
individual cell's own (xi, eta).

No sympy needed anywhere in this script -- RESOLUTIONS can go higher than
the vector-potential counterpart script's sympy-limited range.

Usage: python scripts/line_integral_convergence_scalar_potential.py
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
from generate_cubedsphere_grid import buildCubedSphereGrid, facesAlongPath  # noqa: E402
from scalar_potential_reconstruction_convergence import potential, buildExactEdgeData, fitAlpha  # noqa: E402

PATH_P0 = (-30.0, -20.0)
PATH_P1 = (45.0, 50.0)

RESOLUTIONS = (4, 8, 16, 32, 64, 128)


def computeFlux(M):
    """Build the grid + exact edge data at this M, then PolylineIntegral's flux along PATH_P0->PATH_P1."""
    grid = buildCubedSphereGrid(M)
    data = buildExactEdgeData(grid)

    xyz = numpy.array([(PATH_P0[0], PATH_P0[1], 0.), (PATH_P1[0], PATH_P1[1], 0.)])
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
    pli.computeWeights(xyz, counterclock=False)
    return pli.getIntegral(data, placement=CELL_BY_CELL_DATA)


def main():
    faces = facesAlongPath(PATH_P0, PATH_P1)
    print(f'path {PATH_P0} -> {PATH_P1} crosses cube faces: {sorted(faces)}')
    assert len(faces) >= 2, 'path does not cross multiple cubed-sphere panels'

    exact = potential(PATH_P1[0], PATH_P1[1]) - potential(PATH_P0[0], PATH_P0[1])
    print(f'exact line integral (chi(p1) - chi(p0), no quadrature needed) = {exact:.6f}\n')

    Ms, errs = [], []
    for M in RESOLUTIONS:
        flux = computeFlux(M)
        err = abs(flux - exact)
        Ms.append(M)
        errs.append(err)
        print(f'M={M:4d} (ncells={6 * M * M:7d}): flux={flux:.6f}  err={err:.3e}')

    alpha, prefactor = fitAlpha(Ms, errs)
    print(f'\nfitted error ~ M^-alpha: alpha = {alpha:.3f}')

    plotResult(Ms, errs, alpha, prefactor)
    return alpha


def plotResult(Ms, errs, alpha, prefactor):
    fig, ax = plt.subplots(figsize=(7, 6))
    Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
    ax.loglog(Ms_arr, errs, 'o-', color='tab:blue', linewidth=2,
              label=f'|flux - exact|  (alpha={alpha:.2f})')
    ax.loglog(Ms_arr, prefactor * Ms_arr ** (-alpha), 'k--', linewidth=1,
              label=f'fit: {prefactor:.2e} * M^-{alpha:.2f}')

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel('|line integral error|')
    ax.set_title('Line integral convergence: v = grad(chi)\n'
                  r'$\chi=\sin\lambda\,\cos^2\theta$, path crosses +X -> +Z panels')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'line_integral_convergence_scalar_potential.png'
    fig.savefig(outfile, dpi=150)
    print(f'wrote {outfile}')


if __name__ == '__main__':
    main()
