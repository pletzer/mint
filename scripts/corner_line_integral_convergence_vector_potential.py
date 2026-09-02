#!/usr/bin/env python
"""
Isolates whether the RAW LINE INTEGRAL itself (not divided by branch
length, not compared to a point value) is degraded at a cubed-sphere
3-panel junction -- follow-up to corner_reconstruction_convergence_
vector_potential.py, which found the "+" stencil's reconstructed VECTOR
only converges at alpha~1.29 there (vs alpha~1.9-2.2 at the other 3 corner
types), and hypothesised this is because the mirror-symmetry cancellation
argument (see that script's docstring) breaks down where only 3, not an
even/mirror-paired number of, cells meet.

This script reuses the EXACT SAME "+" stencil branches (same 4 corners,
same halfWidthDeg formula, same TEST_CORNERS/cellIndex) but measures a
DIFFERENT quantity for each branch: |flux - exact_line_integral| (mint's
PolylineIntegral flux vs pole_asymmetric_convergence.exactLineIntegral's
closed-form reference for that SAME short straight path), i.e. exactly the
quantity line_integral_convergence_vector_potential.py measured for its
long, multi-cell path -- now applied to these short, corner-centred,
mostly single-or-few-cell branches instead.

This cleanly separates the two effects bundled together in the corner-
reconstruction script's result:
  (a) the raw line integral's own accuracy along the short path, vs
  (b) "average over a short path ~= point value at its midpoint", an
      independent O(halfWidthDeg^2) = O(h^2) effect from shrinking the
      branch width, unrelated to how well PolylineIntegral itself performs.
If (a) alone is already degraded at the 3-panel junction, that confirms the
corner-reconstruction result traces back to the line integral itself, not
to some artefact of dividing by length or of exactVector's point-value
comparison.

West-east and south-north branches are tracked (and printed) SEPARATELY as
well as combined (Euclidean norm of the two branches' errors, mirroring how
corner_reconstruction_convergence_vector_potential.py combined its u/v
component errors into one vector-norm error) -- the 3-panel junction's
geometry is not symmetric between the two branch directions, so they may
behave differently.

Edge data is EXACT (sympy) as in the corner-reconstruction script; the
reference for each branch is comparatively cheap (ONE exactLineIntegral
call per branch, not per grid edge). Same RESOLUTIONS/runtime constraint:
capped at (4, 8, 16, 32), run in the background.

Usage: python scripts/corner_line_integral_convergence_vector_potential.py
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
from generate_cubedsphere_grid import buildCubedSphereGrid  # noqa: E402
from pole_asymmetric_convergence import exactLineIntegral  # noqa: E402
from vector_potential_reconstruction_convergence import buildExactEdgeData  # noqa: E402
from scalar_potential_reconstruction_convergence import fitAlpha  # noqa: E402
from corner_reconstruction_convergence_vector_potential import (  # noqa: E402
    RESOLUTIONS, TEST_CORNERS, cellIndex, branchIntegral)


def cornerBranchErrors(grid, data, lon0, lat0, halfWidthDeg):
    """
    Raw line-integral error for the west-east and south-north branches
    centred at (lon0, lat0) -- |flux - exact|, NOT divided by length.

    :returns: (err_we, err_sn)
    """
    we_p0, we_p1 = (lon0 - halfWidthDeg, lat0), (lon0 + halfWidthDeg, lat0)
    sn_p0, sn_p1 = (lon0, lat0 - halfWidthDeg), (lon0, lat0 + halfWidthDeg)

    we_flux = branchIntegral(grid, data, we_p0, we_p1)
    we_exact = exactLineIntegral(we_p0, we_p1)
    err_we = abs(we_flux - we_exact)

    sn_flux = branchIntegral(grid, data, sn_p0, sn_p1)
    sn_exact = exactLineIntegral(sn_p0, sn_p1)
    err_sn = abs(sn_flux - sn_exact)

    return err_we, err_sn


def main():
    # errsByLabel[label] = (list of M, list of err_we, list of err_sn, list of combined)
    errsByLabel = {label: ([], [], [], []) for label, _ in TEST_CORNERS}

    for M in RESOLUTIONS:
        print(f'M={M}: building exact (sympy) edge data for {6 * M * M} cells '
              f'({6 * M * M * 4} edges, slow)...', flush=True)
        grid = buildCubedSphereGrid(M)
        data = buildExactEdgeData(grid)
        points = grid.getPoints()

        halfWidthDeg = 0.5 * (90.0 / M)

        for label, cornerFn in TEST_CORNERS:
            faceIdx, I, J, cornerIdx = cornerFn(M)
            lon0, lat0 = points[cellIndex(faceIdx, I, J, M), cornerIdx, :2]

            err_we, err_sn = cornerBranchErrors(grid, data, lon0, lat0, halfWidthDeg)
            combined = numpy.hypot(err_we, err_sn)

            Ms_, wes_, sns_, combs_ = errsByLabel[label]
            Ms_.append(M)
            wes_.append(err_we)
            sns_.append(err_sn)
            combs_.append(combined)

            print(f'  M={M:3d} {label:18s} corner=({lon0:8.3f},{lat0:7.3f})  '
                  f'err_we={err_we:.3e}  err_sn={err_sn:.3e}  combined={combined:.3e}', flush=True)

    print()
    alphas = {}
    for label, (Ms, wes, sns, combs) in errsByLabel.items():
        a_we, _ = fitAlpha(Ms, wes)
        a_sn, _ = fitAlpha(Ms, sns)
        a_comb, _ = fitAlpha(Ms, combs)
        alphas[label] = (a_we, a_sn, a_comb)
        print(f'{label:18s}: alpha_we={a_we:.3f}  alpha_sn={a_sn:.3f}  alpha_combined={a_comb:.3f}')

    plotResult(errsByLabel, alphas)
    return alphas


def plotResult(errsByLabel, alphas):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    colors = plt.cm.tab10.colors

    for i, (label, (Ms, wes, sns, combs)) in enumerate(errsByLabel.items()):
        Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
        a_we, a_sn, a_comb = alphas[label]
        color = colors[i % len(colors)]
        ax.loglog(Ms_arr, combs, 'o-', color=color, linewidth=2,
                   label=f'{label}  (combined alpha={a_comb:.2f})')

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel(r'combined $|$flux $-$ exact$|$ over west-east + south-north branches')
    ax.set_title('Raw line-integral error along short "+" stencil branches\n'
                  r'$\mathbf{u} = \nabla\times(\psi\,\hat{r})$, '
                  r'$\psi=\cos\theta(1+\sin\theta)\cos\lambda$')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'corner_line_integral_convergence_vector_potential.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
