#!/usr/bin/env python
"""
Phase 2 of trajectory_contravariant_PLAN.md: pointwise accuracy of the
CONTRAVARIANT velocity reconstructed inside the cells of an equiangular
cubed sphere -- the quantity a trajectory integrator in parametric space
needs (d xi^i / dt = v^i) -- without any time integration yet.

Arms (all from exact data on the true great-circle cell edges):
  A  Whitney k=1 from edge circulations (covariant v_i), raised with the
     EXACT analytic metric g^ij of the equiangular chart.
  B  Whitney k=1 from the same circulations, raised with the DISCRETE metric
     of the bilinear (xyz) map through the cell's 4 corners -- what a code
     that only knows node coordinates does (mint's getEdgeVectors raises with
     this kind of per-cell bilinear Jacobian). The contravariant components
     are then used as d xi^i/dt in the exact chart.
  C  RT0 / C-grid from edge FLUXES (contravariant density sqrt(g) v^i),
     divided by the exact sqrt(g): no metric inverse needed.
  D  nodal star recovery + bilinear interpolation (nodal_star.py), projected
     onto the exact dual basis.

The error norm is |(v_h^i - v^i) a_i| (physical length of the error in the
contravariant components), normalised by the RMS of |v|. RMS over cells,
split into cell classes (interior / panel-edge / cube-corner cells),
at the centroid and pooled over off-centroid points.

Usage: python contravariant_pointwise_convergence.py [--maxM 128] [--maxMnodal 64]
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
                        whitneyInterp, greatCircleIntegral)
from cs_chart import basis, chartPoint  # noqa: E402
from trajectory_fields import SolidBody, NairLauritzen  # noqa: E402

OFF_CENTRE = ((0.25, 0.25), (0.1, 0.7), (0.9, 0.4), (0.75, 0.05), (0.4, 0.95))
CENTROID = (0.5, 0.5)
ARMS = ('A', 'B', 'C', 'D')
ARM_LABELS = {'A': 'A Whitney + exact g^ij', 'B': 'B Whitney + bilinear g^ij',
              'C': 'C RT0 fluxes / sqrt(g)', 'D': 'D nodal star + bilinear'}
CLASSES = ('interior', 'panel_edge', 'corner')

OMEGA = numpy.array([numpy.sin(0.7) * numpy.cos(0.3),
                     numpy.sin(0.7) * numpy.sin(0.3), numpy.cos(0.7)])
FIELDS = {'solid_body': SolidBody(OMEGA), 'nair_lauritzen': NairLauritzen()}
T_SNAPSHOT = 1.3


def cellClass(M):
    i = numpy.arange(M)
    bi = (i == 0) | (i == M - 1)
    I, J = numpy.meshgrid(bi, bi, indexing='ij')
    return numpy.where(I & J, 'corner', numpy.where(I | J, 'panel_edge', 'interior'))


def panelData(cs, f, fld, t):
    """Edge circulations and fluxes on panel f, oriented along +alpha / +beta."""
    P = cs.nodes[cs.nodeId[f]]                       # (M+1, M+1, 3), [i (alpha), j (beta)]
    vel = (lambda x: fld.velocity(x, t))
    if isinstance(fld, SolidBody):
        circ = fld.edgeIntegral
    else:
        circ = (lambda p, q: greatCircleIntegral(vel, p, q, nGauss=12))
    Ex = circ(P[:-1, :], P[1:, :])                  # (M, M+1): xi-edges (along +alpha)
    Ey = circ(P[:, :-1], P[:, 1:])                  # (M+1, M): eta-edges (along +beta)
    psi = fld.psi(P, t)
    # flux in the +xi direction through eta-edges, in the +eta direction through xi-edges
    Fx = psi[:, :-1] - psi[:, 1:]                   # (M+1, M)
    Fy = psi[1:, :] - psi[:-1, :]                   # (M, M+1)
    return P, Ex, Ey, Fx, Fy


def reconstructPanel(cs, f, data, xi, eta):
    """Arms A, B, C on panel f at cell coordinates (xi, eta). Returns dict of (M, M, 3)."""
    P, Ex, Ey, Fx, Fy = data
    a = cs.alphas
    da = numpy.diff(a)
    AA, BB = numpy.meshgrid(a[:-1] + xi * da, a[:-1] + eta * da, indexing='ij')
    ja, jb = basis(f, AA, BB)
    jx = ja * da[:, None, None]                      # exact d x / d xi
    jy = jb * da[None, :, None]                      # exact d x / d eta
    x = chartPoint(f, AA, BB)
    # covariant cell components (Whitney k=1)
    vxi = (1 - eta) * Ex[:, :-1] + eta * Ex[:, 1:]
    veta = (1 - xi) * Ey[:-1, :] + xi * Ey[1:, :]

    def raiseWith(gx, gy):
        g11 = numpy.sum(gx * gx, -1)
        g12 = numpy.sum(gx * gy, -1)
        g22 = numpy.sum(gy * gy, -1)
        det = g11 * g22 - g12 ** 2
        return (g22 * vxi - g12 * veta) / det, (-g12 * vxi + g11 * veta) / det

    out = {}
    c1, c2 = raiseWith(jx, jy)
    out['A'] = c1[..., None] * jx + c2[..., None] * jy
    # discrete bilinear (xyz, unprojected) Jacobian through the 4 corners
    p00, p10, p01, p11 = P[:-1, :-1], P[1:, :-1], P[:-1, 1:], P[1:, 1:]
    bx = (1 - eta) * (p10 - p00) + eta * (p11 - p01)
    by = (1 - xi) * (p01 - p00) + xi * (p11 - p10)
    c1, c2 = raiseWith(bx, by)
    out['B'] = c1[..., None] * jx + c2[..., None] * jy
    # RT0 / contravariant Piola with exact (signed) Jacobian
    J = numpy.sum(x * numpy.cross(jx, jy), -1)
    fxi = (1 - xi) * Fx[:-1, :] + xi * Fx[1:, :]
    feta = (1 - eta) * Fy[:, :-1] + eta * Fy[:, 1:]
    out['C'] = (fxi / J)[..., None] * jx + (feta / J)[..., None] * jy
    return out, x


def rms(e, mask=None):
    s = numpy.sum(e ** 2, -1)
    if mask is not None:
        s = s[mask]
    return numpy.sqrt(numpy.mean(s))


def fitAlpha(Ms, errs):
    return -numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)[0]


def run(fname, fld, Ms, maxMnodal):
    # errs[(arm, where, cls)] -> list over M ; where in ('centroid', 'off')
    errs = {}
    for M in Ms:
        t0 = time.time()
        cs = CubedSphere(M)
        cls = numpy.stack([cellClass(M)] * 6)
        doNodal = M <= maxMnodal
        if doNodal:
            V = nodalVectors(cs, cs.edgeIntegrals(lambda x: fld.velocity(x, T_SNAPSHOT), nGauss=12))
        datas = [panelData(cs, f, fld, T_SNAPSHOT) for f in range(6)]
        for where, pts in (('centroid', [CENTROID]), ('off', OFF_CENTRE)):
            sq = {arm: [] for arm in ARMS}
            vmag = []
            for (xi, eta) in pts:
                rec = [reconstructPanel(cs, f, datas[f], xi, eta) for f in range(6)]
                x = numpy.stack([r[1] for r in rec])
                exact = fld.velocity(x, T_SNAPSHOT)
                vmag.append(numpy.sum(exact ** 2, -1))
                for arm in 'ABC':
                    sq[arm].append(numpy.stack([r[0][arm] for r in rec]) - exact)
                if doNodal:
                    sq['D'].append(bilinearInterp(cs, V, xi, eta) - exact)
                if M == Ms[0] and isinstance(fld, SolidBody):
                    # cross-check arm A against nodal_star.whitneyInterp
                    wi = whitneyInterp(cs, cs.edgeIntegrals(lambda y: fld.velocity(y)), xi, eta)
                    Aall = numpy.stack([r[0]['A'] for r in rec])
                    assert numpy.allclose(wi, Aall, atol=1e-12), numpy.abs(wi - Aall).max()
            vref = numpy.sqrt(numpy.mean(vmag))
            for arm in ARMS:
                if not sq[arm]:
                    continue
                e = numpy.stack(sq[arm])                       # (npts, 6, M, M, 3)
                for c in CLASSES + ('all',):
                    mask = None if c == 'all' else numpy.broadcast_to(cls == c, e.shape[:-1])
                    errs.setdefault((arm, where, c), []).append(rms(e, mask) / vref)
        print(f'  {fname} M={M:4d} {time.time() - t0:6.1f}s')
    return errs


def report(fname, errs, Ms):
    print(f'\n=== {fname}: relative RMS error of contravariant velocity ===')
    print(f'{"arm":28s} {"where":9s} {"class":11s}' + ''.join(f'M={M:<9d}' for M in Ms) + 'alpha  (last-pair)')
    for arm in ARMS:
        for where in ('centroid', 'off'):
            for c in ('all',) + CLASSES:
                e = errs.get((arm, where, c))
                if not e:
                    continue
                m = Ms[:len(e)]
                loc = numpy.log2(e[-2] / e[-1]) if len(e) > 1 else numpy.nan
                print(f'{ARM_LABELS[arm]:28s} {where:9s} {c:11s}' + ''.join(f'{v:.3e}  ' for v in e)
                      + f'{fitAlpha(m, e):5.2f}  ({loc:.2f})')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--maxM', type=int, default=128)
    parser.add_argument('--maxMnodal', type=int, default=64)
    args = parser.parse_args()
    Ms = [M for M in (8, 16, 32, 64, 128, 256) if M <= args.maxM]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    for row, (fname, fld) in enumerate(FIELDS.items()):
        errs = run(fname, fld, Ms, args.maxMnodal)
        report(fname, errs, Ms)
        for col, where in enumerate(('off', 'centroid')):
            ax = axes[row, col]
            for k, arm in enumerate(ARMS):
                for c, ls in (('interior', '-'), ('corner', ':')):
                    e = errs.get((arm, where, c))
                    if not e:
                        continue
                    m = Ms[:len(e)]
                    ax.loglog(m, e, ls, marker='os^d'[k], color=f'C{k}',
                              label=f'{ARM_LABELS[arm]}, {c} ({fitAlpha(m, e):.2f})')
            ax.set_title(f'{fname}: {"off-centroid points" if where == "off" else "cell centroid"}')
            ax.set_xlabel('M (cells per panel edge)')
            ax.set_ylabel('relative RMS |(v_h^i - v^i) a_i|')
            ax.legend(fontsize=7)
    fig.tight_layout()
    out = HERE / 'contravariant_pointwise_convergence.png'
    fig.savefig(out, dpi=120)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
