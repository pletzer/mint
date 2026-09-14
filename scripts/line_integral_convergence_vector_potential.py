#!/usr/bin/env python
"""
Convergence of mint.PolylineIntegral's line integral of u = curl(psi * r_hat)
along the SAME path as line_integral_convergence_scalar_potential.py
(crosses cube panels, verified via generate_cubedsphere_grid.facesAlongPath),
as the cubed-sphere grid is refined -- the vector-potential counterpart to
that script.

psi(lambda, theta) = cos(theta)*(1+sin(theta))*cos(lambda), the same
pole-asymmetric wavenumber-1 stream function used throughout this codebase
(scripts/pole_asymmetric_convergence.py, sympy-verified there; also
scripts/vector_potential_reconstruction_convergence.py).

Edge data is EXACT here too, per request: buildExactEdgeData (reused
unmodified from vector_potential_reconstruction_convergence.py) integrates
u.dl along every grid edge in closed form via
pole_asymmetric_convergence.exactLineIntegral (sympy, no quadrature) --
unlike the scalar-potential script, this is NOT exact by the fundamental
theorem of calculus (a genuinely rotational field has real circulation, so
its edge integral has to actually be computed, not just differenced from a
potential). "Exact" here specifically means: no quadrature error in
evaluating the integral along the flat-(lon,lat)-chord approximation of
each edge -- it does NOT mean the flat-chord itself matches the true
(closer to great-circle) mesh edge; see that script's module docstring for
the full discussion (empirically confirmed there: quadrature error was
already negligible, so this mainly matters for correctness/rigor, not for
materially changing the numbers).

The REFERENCE value for the WHOLE PATH is comparatively cheap: ONE call to
pole_asymmetric_convergence.exactLineIntegral(PATH_P0, PATH_P1) (same
closed-form sympy method, applied to the path's own endpoints rather than
per grid edge), independent of M.

Cost/scope: exactLineIntegral redoes a fresh sympy integration per edge,
measured at ~60ms/edge for this field (see
vector_potential_reconstruction_convergence.py's module docstring) --
building the edge data for the WHOLE grid at each M (6*M^2*4 edges) is the
dominant cost, same as that script, and is INDEPENDENT of how many
target/path evaluations are done per M (just one line integral here, vs 4
(xi,eta) locations there) -- so RESOLUTIONS is capped at (4, 8, 16, 32),
expected to take on the order of ~10-35 minutes total; run this in the
background.

Usage: python scripts/line_integral_convergence_vector_potential.py
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
from pole_asymmetric_convergence import exactLineIntegral  # noqa: E402
from vector_potential_reconstruction_convergence import buildExactEdgeData  # noqa: E402
from scalar_potential_reconstruction_convergence import fitAlpha  # noqa: E402 (generic, field-agnostic helper)

PATH_P0 = (-30.0, -20.0)
PATH_P1 = (45.0, 50.0)

RESOLUTIONS = (4, 8, 16, 32)  # not higher: exactLineIntegral is per-edge sympy, see module docstring


def computeFlux(M):
    """Build the grid + exact (sympy) edge data at this M, then PolylineIntegral's flux along the path."""
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

    exact = exactLineIntegral(PATH_P0, PATH_P1)
    print(f'exact line integral (sympy, flat (lon,lat) chord) = {exact:.6f}\n')

    Ms, errs = [], []
    for M in RESOLUTIONS:
        print(f'M={M}: building exact (sympy) edge data for {6 * M * M} cells '
              f'({6 * M * M * 4} edges, slow)...', flush=True)
        flux = computeFlux(M)
        err = abs(flux - exact)
        Ms.append(M)
        errs.append(err)
        print(f'M={M:4d} (ncells={6 * M * M:7d}): flux={flux:.6f}  err={err:.3e}', flush=True)

    alpha, prefactor = fitAlpha(Ms, errs)
    print(f'\nfitted error ~ M^-alpha: alpha = {alpha:.3f}')

    plotResult(Ms, errs, alpha, prefactor)
    return alpha


def plotResult(Ms, errs, alpha, prefactor):
    fig, ax = plt.subplots(figsize=(7, 6))
    Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
    ax.loglog(Ms_arr, errs, 'o-', color='tab:red', linewidth=2,
              label=f'|flux - exact|  (alpha={alpha:.2f})')
    ax.loglog(Ms_arr, prefactor * Ms_arr ** (-alpha), 'k--', linewidth=1,
              label=f'fit: {prefactor:.2e} * M^-{alpha:.2f}')

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel('|line integral error|')
    ax.set_title('Line integral convergence: ' + r'$\mathbf{u} = \nabla\times(\psi\,\hat{r})$'
                  + '\n' + r'$\psi=\cos\theta(1+\sin\theta)\cos\lambda$, path crosses +X -> +Z panels')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'line_integral_convergence_vector_potential.png'
    fig.savefig(outfile, dpi=150)
    print(f'wrote {outfile}')


if __name__ == '__main__':
    main()
