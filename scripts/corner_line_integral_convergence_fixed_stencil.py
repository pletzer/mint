#!/usr/bin/env python
"""
Second version of corner_line_integral_convergence_vector_potential.py,
fixing two issues raised in the conversation this came out of:

1. FIXED (not M-dependent) stencil size. The original script's branches
   shrank with the grid (halfWidthDeg = 0.5*(90/M)), so the "exact" line
   integral being estimated ALSO shrank like O(h) as M grew -- meaning its
   headline "alpha" numbers were fit on ABSOLUTE error while the quantity
   itself was vanishing at a comparable rate, conflating "the signal is
   getting smaller" with "the error is getting smaller". A genuinely
   flat-in-M relative error near the 3-panel junction's west-east branch
   was misread as a middling alpha~1.08 there. HALF_WIDTH_DEG below is
   fixed across all M -- the same branch, refined underneath by an
   ever-finer grid, exactly mirroring line_integral_convergence_vector_
   potential.py's long, fixed panel-crossing path, just centred on these 4
   special corner locations instead. With a fixed target, ABSOLUTE error is
   directly meaningful again (no shrinking-signal confound); relative error
   is also tracked/printed for comparison.

2. The vector potential's own longitude origin is shifted by LON_SHIFT_DEG
   (a "deliberately non-round angle", matching generate_cubedsphere_grid.py's
   own rotationDeg=23.7 precedent for the same reason). uDotDl's dtheta-
   coefficient (the south-north branch's sensitivity) carries a sin(lambda)
   factor that is EXACTLY zero at lambda=0 and lambda=180 -- precisely the
   'interior' and 'near dateline' test corners -- which is why the previous
   script's south-north errors there were ~1e-16: not a genuine numerical
   cancellation, just psi being locally insensitive to theta at those two
   longitudes, an accident of the corners' own (grid-motivated) round-number
   choice aligning with the field's (unrelated) zero-crossing. Shifting the
   field's longitude origin by a generic, non-round amount decouples the two,
   so no branch is masked by an unrelated symmetry.

Both fixes are independent of the RESOLUTIONS/RankWarning issue: shifting
does NOT change which corners are tested (still 'interior', '2-panel edge',
'3-panel junction', 'near dateline', reused directly from
corner_reconstruction_convergence_vector_potential.TEST_CORNERS) or where
they sit -- only which phase of psi is evaluated there, and edge data is
rebuilt (uDotDlShifted/exactLineIntegralShifted/buildExactEdgeDataShifted
below) to stay self-consistent with that shifted field. Same cost/
RESOLUTIONS constraint as before: exact (sympy) edge data, capped at
(4, 8, 16, 32), run in the background (~10-35 minutes).

Usage: python scripts/corner_line_integral_convergence_fixed_stencil.py
"""
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import NUM_EDGES_PER_QUAD

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from generate_cubedsphere_grid import buildCubedSphereGrid  # noqa: E402
from pole_asymmetric_convergence import uDotDl, exactLineIntegral  # noqa: E402
from scalar_potential_reconstruction_convergence import fitAlpha  # noqa: E402
from corner_reconstruction_convergence_vector_potential import (  # noqa: E402
    RESOLUTIONS, TEST_CORNERS, cellIndex, branchIntegral)

DEG2RAD = numpy.pi / 180.0

LON_SHIFT_DEG = 23.7  # non-round, see module docstring's point (2)
HALF_WIDTH_DEG = 5.0  # FIXED across all M, see module docstring's point (1)


# ---------------------------------------------------------------------------
# The longitude-shifted vector potential: same psi = cos(theta)(1+sin(theta))
# cos(lambda - LON_SHIFT_DEG), just with its own phase offset. Re-derived
# here (rather than mutating the shared, unshifted helpers other scripts
# rely on) so those scripts' own behaviour/output stays untouched.
# ---------------------------------------------------------------------------
def exactLineIntegralShifted(p0, p1):
    """exactLineIntegral with both endpoints' longitude shifted by LON_SHIFT_DEG."""
    return exactLineIntegral((p0[0] - LON_SHIFT_DEG, p0[1]), (p1[0] - LON_SHIFT_DEG, p1[1]))


def buildExactEdgeDataShifted(grid):
    """Like vector_potential_reconstruction_convergence.buildExactEdgeData, for the shifted field."""
    ncells = grid.getNumberOfCells()
    points = grid.getPoints()
    data = numpy.zeros((ncells, NUM_EDGES_PER_QUAD))
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        for icell in range(ncells):
            data[icell, i0] = sign * exactLineIntegralShifted(
                points[icell, i0, :2], points[icell, i1, :2])
    return data


def cornerBranchErrors(grid, data, lon0, lat0, halfWidthDeg):
    """
    Raw line-integral error (absolute AND relative) for the west-east and
    south-north branches centred at (lon0, lat0), FIXED width, shifted field.

    :returns: (err_we, rel_we, err_sn, rel_sn)
    """
    we_p0, we_p1 = (lon0 - halfWidthDeg, lat0), (lon0 + halfWidthDeg, lat0)
    sn_p0, sn_p1 = (lon0, lat0 - halfWidthDeg), (lon0, lat0 + halfWidthDeg)

    we_flux = branchIntegral(grid, data, we_p0, we_p1)
    we_exact = exactLineIntegralShifted(we_p0, we_p1)
    err_we = abs(we_flux - we_exact)
    rel_we = err_we / abs(we_exact) if abs(we_exact) > 1.e-12 else float('nan')

    sn_flux = branchIntegral(grid, data, sn_p0, sn_p1)
    sn_exact = exactLineIntegralShifted(sn_p0, sn_p1)
    err_sn = abs(sn_flux - sn_exact)
    rel_sn = err_sn / abs(sn_exact) if abs(sn_exact) > 1.e-12 else float('nan')

    return err_we, rel_we, err_sn, rel_sn


def main():
    # errsByLabel[label] = (Ms, err_we, rel_we, err_sn, rel_sn, combined_abs, combined_rel)
    errsByLabel = {label: ([], [], [], [], [], [], []) for label, _ in TEST_CORNERS}

    for M in RESOLUTIONS:
        print(f'M={M}: building exact (sympy) edge data for {6 * M * M} cells '
              f'({6 * M * M * 4} edges, slow)...', flush=True)
        grid = buildCubedSphereGrid(M)
        data = buildExactEdgeDataShifted(grid)
        points = grid.getPoints()

        for label, cornerFn in TEST_CORNERS:
            faceIdx, I, J, cornerIdx = cornerFn(M)
            lon0, lat0 = points[cellIndex(faceIdx, I, J, M), cornerIdx, :2]

            err_we, rel_we, err_sn, rel_sn = cornerBranchErrors(grid, data, lon0, lat0, HALF_WIDTH_DEG)
            combined_abs = numpy.hypot(err_we, err_sn)
            combined_rel = numpy.hypot(rel_we, rel_sn)

            Ms_, wes_, rwes_, sns_, rsns_, cabs_, crel_ = errsByLabel[label]
            Ms_.append(M)
            wes_.append(err_we)
            rwes_.append(rel_we)
            sns_.append(err_sn)
            rsns_.append(rel_sn)
            cabs_.append(combined_abs)
            crel_.append(combined_rel)

            print(f'  M={M:3d} {label:18s} corner=({lon0:8.3f},{lat0:7.3f})  '
                  f'err_we={err_we:.3e} (rel={rel_we:.3e})  '
                  f'err_sn={err_sn:.3e} (rel={rel_sn:.3e})  '
                  f'combined_abs={combined_abs:.3e}', flush=True)

    print()
    alphas = {}
    for label, (Ms, wes, rwes, sns, rsns, cabs, crel) in errsByLabel.items():
        a_we, _ = fitAlpha(Ms, wes)
        a_rel_we, _ = fitAlpha(Ms, rwes)
        a_sn, _ = fitAlpha(Ms, sns)
        a_rel_sn, _ = fitAlpha(Ms, rsns)
        a_comb, _ = fitAlpha(Ms, cabs)
        alphas[label] = dict(we=a_we, rel_we=a_rel_we, sn=a_sn, rel_sn=a_rel_sn, combined=a_comb)
        print(f'{label:18s}: alpha_we={a_we:.3f} (rel={a_rel_we:.3f})  '
              f'alpha_sn={a_sn:.3f} (rel={a_rel_sn:.3f})  alpha_combined={a_comb:.3f}')

    plotResult(errsByLabel, alphas)
    return alphas


def plotResult(errsByLabel, alphas):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    colors = plt.cm.tab10.colors

    for i, (label, (Ms, wes, rwes, sns, rsns, cabs, crel)) in enumerate(errsByLabel.items()):
        Ms_arr = numpy.asarray(Ms, dtype=numpy.float64)
        a_comb = alphas[label]['combined']
        color = colors[i % len(colors)]
        ax.loglog(Ms_arr, cabs, 'o-', color=color, linewidth=2,
                   label=f'{label}  (combined alpha={a_comb:.2f})')

    ax.set_xlabel('M  (cells per cubed-sphere panel edge; 6*M^2 cells total)')
    ax.set_ylabel(r'combined $|$flux $-$ exact$|$, FIXED-width branches')
    ax.set_title(f'Line-integral error, fixed {2*HALF_WIDTH_DEG:.0f} deg branches, '
                  f'longitude-shifted field\n'
                  r'$\mathbf{u} = \nabla\times(\psi\,\hat{r})$, '
                  r'$\psi=\cos\theta(1+\sin\theta)\cos(\lambda-23.7^\circ)$')
    ax.legend(fontsize='small')
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()

    outfile = HERE / 'corner_line_integral_convergence_fixed_stencil.png'
    fig.savefig(outfile, dpi=150)
    print(f'\nwrote {outfile}')


if __name__ == '__main__':
    main()
