#!/usr/bin/env python
"""
Phase 4 of trajectory_contravariant_PLAN.md: convergence of trajectories
integrated cell by cell in parametric space (trajectory_integrator.py)
with the contravariant velocity from reconstruction arms A-E, on an
equiangular cubed sphere with M = 8 ... 128 cells per panel edge.

Cases
  sb_generic   solid-body rotation about an axis tilted 0.7 rad from z,
               random seeds, one revolution (T = 2 pi). Error: endpoint
               and max over 8 intermediate times vs the exact rotation.
  sb_aligned25 rotation about z, seeds near the equator at latitude
  sb_aligned50 q * (pi/2)/M (q = 0.25, 0.50), i.e. at a FIXED sub-cell offset
               from the equator grid line, so paths run (nearly) along
               grid lines on the equatorial panels: tests whether the
               O(h) Whitney error still cancels when it is not averaged
               across cells.
  nl           Nair-Lauritzen deformational flow, random seeds, error of
               the return to the start at t = T = 5.

Error = great-circle distance (radians on the unit sphere); RMS and max
over seeds, plus median and 99th percentile (Whitney arms can trap a few
particles on edges where the normal velocity converges from both sides --
a Filippov sliding mode -- which shows up in max/RMS but not the median).
For sb_generic and nl the RMS is also split by whether the
path ever entered a cube-corner cell. Time step dt = CFL * (pi/2)/M / vmax,
arm E (exact velocity) measures the time-integration + hand-over error.

Usage: python trajectory_convergence_cubedsphere.py [--maxM 128] [--nseeds 400]
                                                    [--cases sb_generic,nl] [--arms ABCDWRSE] [--tag _nl256]
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
from nodal_star import arcLength  # noqa: E402
from trajectory_fields import SolidBody, NairLauritzen  # noqa: E402
from trajectory_integrator import Tracer  # noqa: E402

ARM_LABELS = {'A': 'A Whitney + exact g^ij', 'B': 'B Whitney + bilinear g^ij',
              'C': 'C RT0 fluxes / sqrt(g)', 'D': 'D nodal star + bilinear',
              'W': 'W Whitney k=2 + exact g^ij',
              'R': 'R mint k=2 2x2 + exact g^ij',
              'S': 'S mint k=2 2x2 + biq. g^ij',
              'E': 'E exact v (time error)'}
OMEGA_TILT = numpy.array([numpy.sin(0.7) * numpy.cos(0.3),
                          numpy.sin(0.7) * numpy.sin(0.3), numpy.cos(0.7)])
CFL = 0.5
NOUT = 8


def randomSeeds(n, seed=3):
    x = numpy.random.default_rng(seed).normal(size=(n, 3))
    return x / numpy.linalg.norm(x, axis=-1, keepdims=True)


def equatorSeeds(n, M, q, seed=4):
    lon = numpy.random.default_rng(seed).uniform(0, 2 * numpy.pi, n)
    lat = q * (numpy.pi / 2) / M
    return numpy.stack([numpy.cos(lat) * numpy.cos(lon), numpy.cos(lat) * numpy.sin(lon),
                        numpy.full(n, numpy.sin(lat))], -1)


def caseSetup(case, M, nseeds):
    """Returns field, seeds, final time, output times, exact(x0, t) or None, vmax."""
    if case == 'sb_generic':
        f = SolidBody(OMEGA_TILT)
        T = 2 * numpy.pi
        return f, randomSeeds(nseeds), T, f.trajectory, 1.
    if case.startswith('sb_aligned'):
        f = SolidBody([0., 0., 1.])
        q = int(case[len('sb_aligned'):]) / 100.
        return f, equatorSeeds(nseeds, M, q), 2 * numpy.pi, f.trajectory, 1.
    if case == 'nl':
        f = NairLauritzen()
        return f, randomSeeds(nseeds), f.T, None, f.k + 2 * numpy.pi / f.T
    raise ValueError(case)


def runCase(case, Ms, nseeds, arms):
    res = {}       # (arm, metric) -> list over M
    for M in Ms:
        fld, x0, T, exact, vmax = caseSetup(case, M, nseeds)
        t0 = time.time()
        tr = Tracer(M, fld, arms=arms)
        tbuild = time.time() - t0
        dt = CFL * (numpy.pi / 2) / M / vmax
        line = f'  {case:13s} M={M:4d} build {tbuild:5.1f}s'
        for arm in arms:
            t0 = time.time()
            tr.nSlide = 0
            st = tr.locate(x0.copy())
            visited = numpy.zeros(len(x0), bool)
            maxErr = numpy.zeros(len(x0))
            touts = numpy.linspace(0, T, NOUT + 1)[1:]
            tprev = 0.
            for tk in touts:
                tr.advance(arm, st, tprev, tk, dt, visited)
                tprev = tk
                if exact is not None:
                    maxErr = numpy.maximum(maxErr, arcLength(tr.position(st), exact(x0, tk)))
            ref = exact(x0, T) if exact is not None else x0
            err = arcLength(tr.position(st), ref)
            if exact is None:
                maxErr = err
            for key, val in (('rms', numpy.sqrt(numpy.mean(err ** 2))), ('max', err.max()),
                             ('median', numpy.median(err)), ('p99', numpy.percentile(err, 99)),
                             ('rms_path_max', numpy.sqrt(numpy.mean(maxErr ** 2))),
                             ('rms_corner', numpy.sqrt(numpy.mean(err[visited] ** 2)) if visited.any() else numpy.nan),
                             ('rms_nocorner', numpy.sqrt(numpy.mean(err[~visited] ** 2)) if (~visited).any() else numpy.nan),
                             ('frac_corner', visited.mean()), ('slides', tr.nSlide)):
                res.setdefault((arm, key), []).append(val)
            line += f' | {arm} {time.time() - t0:5.1f}s'
        print(line, flush=True)
    return res


def fitAlpha(Ms, e):
    e = numpy.asarray(e, float)
    ok = numpy.isfinite(e) & (e > 0)
    if ok.sum() < 2:
        return numpy.nan
    return -numpy.polyfit(numpy.log(numpy.asarray(Ms)[ok]), numpy.log(e[ok]), 1)[0]


def report(case, res, Ms, arms):
    print(f'\n=== {case}: great-circle position error (radians) ===')
    print(f'{"arm":27s} {"metric":13s}' + ''.join(f'M={M:<9d}' for M in Ms) + 'alpha (last pair)')
    for arm in arms:
        for key in ('rms', 'median', 'p99', 'max', 'rms_path_max', 'rms_nocorner', 'rms_corner'):
            e = res[(arm, key)]
            loc = numpy.log2(e[-2] / e[-1]) if len(e) > 1 else numpy.nan
            print(f'{ARM_LABELS[arm]:27s} {key:13s}' + ''.join(f'{v:.3e}  ' for v in e)
                  + f'{fitAlpha(Ms, e):5.2f} ({loc:.2f})')
    print('fraction of paths entering a cube-corner cell: '
          + ' '.join(f'{v:.2f}' for v in res[(arms[0], 'frac_corner')]))
    print('slide hand-overs per arm: ' + '; '.join(
        f'{arm}: ' + ','.join(str(int(v)) for v in res[(arm, 'slides')]) for arm in arms))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--maxM', type=int, default=128)
    parser.add_argument('--nseeds', type=int, default=400)
    parser.add_argument('--cases', default='sb_generic,sb_aligned25,sb_aligned50,nl')
    parser.add_argument('--arms', default='ABCDE')
    parser.add_argument('--tag', default='', help='suffix for the output png')
    args = parser.parse_args()
    Ms = [M for M in (8, 16, 32, 64, 128, 256) if M <= args.maxM]
    cases = args.cases.split(',')
    arms = tuple(args.arms)

    fig, axes = plt.subplots(1, len(cases), figsize=(5.5 * len(cases), 5), squeeze=False)
    for ax, case in zip(axes[0], cases):
        res = runCase(case, Ms, args.nseeds, arms)
        report(case, res, Ms, arms)
        for k, arm in enumerate(arms):
            e = res[(arm, 'rms')]
            ax.loglog(Ms, e, marker='os^dpv<>'[k], color=f'C{k}',
                      label=f'{ARM_LABELS[arm]} ({fitAlpha(Ms, e):.2f})')
        ax.loglog(Ms, res[(arms[0], 'rms')][0] * (Ms[0] / numpy.array(Ms, float)) ** 2, 'k:', lw=0.8,
                  label='M^-2')
        ax.set_title(case)
        ax.set_xlabel('M (cells per panel edge)')
        ax.set_ylabel('RMS endpoint error (rad)')
        ax.legend(fontsize=7)
    fig.tight_layout()
    out = HERE / f'trajectory_convergence_cubedsphere{args.tag}.png'
    fig.savefig(out, dpi=120)
    print(f'\nwrote {out}')


if __name__ == '__main__':
    main()
