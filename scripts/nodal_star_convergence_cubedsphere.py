#!/usr/bin/env python
"""
Convergence of nodal star-stencil recovery + bilinear interpolation
(scripts/nodal_star.py) versus lowest-order Whitney reconstruction, on an
equiangular cubed sphere, from the SAME exact great-circle edge integrals.

Reports, for M = RESOLUTIONS cells per panel edge:
  * nodal errors (max) per node class: interior, panel edge, cube corner,
    plus cube corners with the O(h) single-segment 3-star (rayOrder=1);
  * RMS error over all cells at several (xi, eta) for nodal+bilinear and
    for Whitney (whose 2nd order should only show up at the centroid).
Fits error ~ M^-alpha and writes nodal_star_convergence_cubedsphere.png.

--jitter J randomly (antisymmetrically) displaces the equiangular grid
coordinates by up to J times the spacing, so that neighbouring edges have
O(1) length ratios at every M -- a stress test of the unequal-arm weights
(on the plain equiangular grid neighbouring edges differ only by O(h)).

Usage: python scripts/nodal_star_convergence_cubedsphere.py [--jitter 0.4]
"""
import argparse
import sys
import time
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))
from nodal_star import (CubedSphere, nodalVectors, bilinearInterp,  # noqa: E402
                        whitneyInterp, cellTargets, testField)

RESOLUTIONS = (8, 16, 32, 64, 128)
XI_ETA_POINTS = ((0.5, 0.5), (0.25, 0.25), (0.1, 0.7), (0.9, 0.4))


def fitAlpha(Ms, errs):
    return -numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)[0]


def rms(e):
    return numpy.sqrt(numpy.mean(numpy.sum(e ** 2, -1)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--jitter', type=float, default=0.0)
    args = parser.parse_args()
    nodal = {k: [] for k in ('interior', 'panel_edge', 'corner', 'corner_rayOrder1')}
    interp = {('nodal', p): [] for p in XI_ETA_POINTS}
    interp.update({('whitney', p): [] for p in XI_ETA_POINTS})

    for M in RESOLUTIONS:
        t0 = time.time()
        cs = CubedSphere(M, jitter=args.jitter)
        assert len(cs.nodes) == 6 * M * M + 2, len(cs.nodes)
        assert numpy.sum(cs.nodeClass == 'corner') == 8
        I = cs.edgeIntegrals(testField)
        V = nodalVectors(cs, I)
        err = numpy.linalg.norm(V - testField(cs.nodes), axis=-1)
        for k in ('interior', 'panel_edge', 'corner'):
            nodal[k].append(err[cs.nodeClass == k].max())
        corners = numpy.where(cs.nodeClass == 'corner')[0]
        V1 = nodalVectors(cs, I, rayOrder=1, nodes=corners)
        nodal['corner_rayOrder1'].append(
            numpy.linalg.norm(V1[corners] - testField(cs.nodes[corners]), axis=-1).max())

        for p in XI_ETA_POINTS:
            x, _, _ = cellTargets(cs, *p)
            exact = testField(x)
            interp[('nodal', p)].append(rms(bilinearInterp(cs, V, *p) - exact))
            interp[('whitney', p)].append(rms(whitneyInterp(cs, I, *p) - exact))
        print(f'M={M:4d} done in {time.time() - t0:.1f}s')

    Ms = numpy.array(RESOLUTIONS, float)
    print('\nnodal max error                  ' + ' '.join(f'M={M:<9d}' for M in RESOLUTIONS) + ' alpha')
    for k, e in nodal.items():
        print(f'  {k:30s} ' + ' '.join(f'{v:.3e}  ' for v in e) + f'{fitAlpha(Ms, e):.2f}')
    print('\ncell RMS error at (xi, eta)')
    for (meth, p), e in interp.items():
        print(f'  {meth:8s} {str(p):12s}         ' + ' '.join(f'{v:.3e}  ' for v in e)
              + f'{fitAlpha(Ms, e):.2f}')

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for k, e in nodal.items():
        ax.loglog(Ms, e, 'o-', label=f'{k} (alpha={fitAlpha(Ms, e):.2f})')
    ax.set_title('Nodal recovery: max error per node class')
    ax.set_xlabel('M (cells per panel edge)')
    ax.set_ylabel('max |v - v_exact|')
    ax.legend(fontsize=8)
    ax = axes[1]
    for (meth, p), e in interp.items():
        ax.loglog(Ms, e, ('o-' if meth == 'nodal' else 's--'),
                  label=f'{meth} {p} (alpha={fitAlpha(Ms, e):.2f})')
    ax.set_title('Interpolation inside cells: RMS error')
    ax.set_xlabel('M (cells per panel edge)')
    ax.set_ylabel('RMS |v - v_exact|')
    ax.legend(fontsize=7)
    fig.tight_layout()
    suffix = f'_jitter{args.jitter:g}' if args.jitter > 0 else ''
    out = HERE / f'nodal_star_convergence_cubedsphere{suffix}.png'
    fig.savefig(out, dpi=120)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
