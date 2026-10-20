#!/usr/bin/env python
"""
Visualise trajectories integrated by trajectory_integrator.py with the
different contravariant reconstructions, against the exact-velocity path.

Seeds: NRANDOM random seeds plus the NWORST seeds with the largest Whitney
(arm A) return error among 400 random seeds, so the problematic paths are
always included.

Writes
  trajectories_<case>_M<M>.png   lon-lat overview + zoom on the worst path
  trajectories_<case>_M<M>.html  interactive 3D (plotly): sphere, cubed-sphere
                                 cells, one toggleable trace group per arm

Usage: python plot_trajectories.py [--case nl|sb_generic] [--M 16] [--arms ACDWRE]
"""
import argparse
import sys
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))
from nodal_star import arcLength, facePoint  # noqa: E402
from trajectory_integrator import Tracer  # noqa: E402
from trajectory_convergence_cubedsphere import caseSetup, randomSeeds, CFL, ARM_LABELS  # noqa: E402

NRANDOM = 8
NWORST = 4
NOUT = 240
COLORS = {'A': '#1f77b4', 'B': '#9467bd', 'C': '#ff7f0e', 'D': '#2ca02c', 'W': '#d62728', 'R': '#8c564b', 'E': 'black'}


def integratePaths(tr, arm, x0, T, dt):
    st = tr.locate(x0.copy())
    path = [x0.copy()]
    tprev = 0.
    for tk in numpy.linspace(0, T, NOUT + 1)[1:]:
        tr.advance(arm, st, tprev, tk, dt)
        tprev = tk
        path.append(tr.position(st))
    return numpy.array(path)            # (NOUT+1, nseeds, 3)


def lonlat(x):
    return (numpy.degrees(numpy.arctan2(x[..., 1], x[..., 0])),
            numpy.degrees(numpy.arcsin(numpy.clip(x[..., 2], -1, 1))))


def plotLonLat(ax, x, **kw):
    lon, lat = lonlat(x)
    lon = lon.astype(float).copy()
    lon[1:][numpy.abs(numpy.diff(lon)) > 180] = numpy.nan   # break at the dateline
    return ax.plot(lon, lat, **kw)


def gridLines(M, nsub=8):
    """Cubed-sphere cell edges as a list of (n, 3) polylines, and the panel edges."""
    a = numpy.linspace(-numpy.pi / 4, numpy.pi / 4, M + 1)
    s = numpy.linspace(-numpy.pi / 4, numpy.pi / 4, M * nsub + 1)
    cells, panels = [], []
    for f in range(6):
        for k, ak in enumerate(a):
            for line in (facePoint(f, numpy.full_like(s, ak), s), facePoint(f, s, numpy.full_like(s, ak))):
                (panels if k in (0, M) else cells).append(line)
    return cells, panels


def zoomPanel(ax, paths, arms, x0, k, cells, panels, exact, endZoom=False, legend=True):
    """
    Orthographic view of path k centred on the path (or, with endZoom or for
    a non-localised path such as a full great-circle orbit, on the start
    point, framing the start and all end points).
    """
    c = paths['E'][:, k].mean(0) if 'E' in paths else paths[arms[0]][:, k].mean(0)
    localised = numpy.linalg.norm(c) > 0.6 and not endZoom
    if not localised:
        c = x0[k].copy()
    c = c / numpy.linalg.norm(c)
    e1 = numpy.cross([0., 0., 1.], c) if abs(c[2]) < 0.99 else numpy.array([1., 0., 0.])
    e1 /= numpy.linalg.norm(e1)
    e2 = numpy.cross(c, e1)

    def ortho(x, **kw):
        x = numpy.asarray(x, float).reshape(-1, 3)
        u, v = numpy.degrees(x @ e1), numpy.degrees(x @ e2)
        hidden = (x @ c) < 0.2
        u[hidden] = numpy.nan
        v[hidden] = numpy.nan
        return ax.plot(u, v, **kw)
    for ln in cells:
        ortho(ln, color='0.85', lw=0.5)
    for ln in panels:
        ortho(ln, color='0.45', lw=1.0)
    for arm in arms:
        ortho(paths[arm][:, k], color=COLORS[arm], lw=1.8 if arm == 'E' else 1.2,
              ls='--' if arm == 'E' else '-', marker='.', ms=3, label=ARM_LABELS[arm], zorder=3 if arm == 'E' else 2)
        ortho(paths[arm][-1, k], marker='s', color=COLORS[arm], ms=8, ls='none', zorder=4)
    ortho(x0[k], marker='*', color='r', ms=15, ls='none', zorder=5,
          label='start (= exact end)' if exact is None else 'start')
    if localised:
        allp = numpy.concatenate([paths[a][:, k] for a in arms])
        pad = 4
    else:
        # frame the start and the end points of the O(h^2) arms (D, W, R) only;
        # k=1 / RT0 end points may fall outside -- their squares then sit off-panel
        frame = [a for a in arms if a in 'DWRE'] or list(arms)
        allp = numpy.stack([x0[k]] + [paths[a][-1, k] for a in frame])
    U, V = numpy.degrees(allp @ e1), numpy.degrees(allp @ e2)
    if localised:
        ax.set_xlim(U.min() - pad, U.max() + pad)
        ax.set_ylim(V.min() - pad, V.max() + pad)
    else:
        half = 0.5 * max(numpy.ptp(U), numpy.ptp(V), 1.0) * 1.6
        uc, vc = 0.5 * (U.min() + U.max()), 0.5 * (V.min() + V.max())
        ax.set_xlim(uc - half, uc + half)
        ax.set_ylim(vc - half, vc + half)
    ax.set_aspect('equal')
    lo, la = lonlat(c)
    what = 'path' if localised else 'start/end region (framed on D, W ends)'
    ax.set_title(f'{what}, centred at ({lo:.0f}E, {la:.0f}N)\n'
                 f'dots every T/{NOUT}, squares = end', fontsize=9)
    ax.set_xlabel('degrees (projected)')
    if legend:
        ax.legend(fontsize=7)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', default='nl')
    parser.add_argument('--M', type=int, default=16)
    parser.add_argument('--arms', default='ACDWE')
    args = parser.parse_args()
    M, arms = args.M, tuple(args.arms)

    fld, _, T, exact, vmax = caseSetup(args.case, M, 1)
    dt = CFL * (numpy.pi / 2) / M / vmax
    tr = Tracer(M, fld, arms=arms)

    # pick the worst Whitney paths among 400 seeds
    cand = randomSeeds(400, seed=11)
    st = tr.locate(cand.copy())
    tr.advance('A', st, 0., T, dt)
    ref = exact(cand, T) if exact is not None else cand
    errA = arcLength(tr.position(st), ref)
    worst = numpy.argsort(errA)[::-1][:NWORST]
    x0 = numpy.concatenate([cand[worst], randomSeeds(NRANDOM, seed=12)])
    print('worst Whitney (A) return errors (rad):', ' '.join(f'{e:.3f}' for e in errA[worst]))

    paths, slides = {}, {}
    for arm in arms:
        tr.nSlide = 0
        paths[arm] = integratePaths(tr, arm, x0, T, dt)
        slides[arm] = tr.nSlide
        ref = exact(x0, T) if exact is not None else x0
        e = arcLength(paths[arm][-1], ref)
        print(f'{ARM_LABELS[arm]:27s} return error: worst-{NWORST} '
              + ' '.join(f'{v:.3f}' for v in e[:NWORST]) + f' | random rms {numpy.sqrt(numpy.mean(e[NWORST:] ** 2)):.2e}'
              + f' | slides {slides[arm]}')

    cells, panels = gridLines(M, nsub=4)
    tag = f'{args.case}_M{M}'

    # ---------------- static lon-lat figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={'width_ratios': [1.6, 1]})
    for ax in axes:
        for ln in panels:
            plotLonLat(ax, ln, color='0.55', lw=0.8)
    for ln in cells:
        plotLonLat(axes[1], ln, color='0.85', lw=0.5)
    for arm in arms:
        for k in range(len(x0)):
            plotLonLat(axes[0], paths[arm][:, k], color=COLORS[arm], lw=1.4 if arm == 'E' else 1.0,
                       ls='--' if arm == 'E' else '-', label=ARM_LABELS[arm] if k == 0 else None)
    lon0, lat0 = lonlat(x0)
    axes[0].plot(lon0[:NWORST], lat0[:NWORST], 'r*', ms=12, label=f'start, {NWORST} worst Whitney paths')
    axes[0].plot(lon0[NWORST:], lat0[NWORST:], 'ko', ms=5, label='start, random seeds')
    axes[0].set_xlim(-180, 180)
    axes[0].set_ylim(-90, 90)
    axes[0].set_xlabel('longitude')
    axes[0].set_ylabel('latitude')
    axes[0].set_title(f'{args.case}, M={M}: all paths (exact v dashed)')
    axes[0].legend(fontsize=8, loc='lower left')

    # zoom on the worst path: orthographic projection centred on that path
    zoomPanel(axes[1], paths, arms, x0, 0, cells, panels, exact)
    fig.tight_layout()
    out = HERE / f'trajectories_{tag}.png'
    fig.savefig(out, dpi=130)
    print(f'wrote {out}')

    # ---------------- worst paths: full path (top) and start/end region (bottom) ----------------
    fig, axes = plt.subplots(2, NWORST, figsize=(5 * NWORST, 10))
    for k in range(NWORST):
        zoomPanel(axes[0, k], paths, arms, x0, k, cells, panels, exact, legend=(k == 0))
        zoomPanel(axes[1, k], paths, arms, x0, k, cells, panels, exact, endZoom=True, legend=False)
    fig.suptitle(f'{args.case}, M={M}: the {NWORST} paths with the largest Whitney k=1 return error '
                 f'(top: whole path; bottom: start/end region)')
    fig.tight_layout()
    out = HERE / f'trajectories_{tag}_worst.png'
    fig.savefig(out, dpi=110)
    print(f'wrote {out}')

    # ---------------- interactive 3D ----------------
    import plotly.graph_objects as go
    traces = []
    u, v = numpy.mgrid[0:2 * numpy.pi:80j, 0:numpy.pi:40j]
    r = 0.995
    traces.append(go.Surface(x=r * numpy.cos(u) * numpy.sin(v), y=r * numpy.sin(u) * numpy.sin(v),
                             z=r * numpy.cos(v), showscale=False, opacity=1.0, hoverinfo='skip',
                             colorscale=[[0, '#f4f4f4'], [1, '#f4f4f4']], name='sphere'))

    def polyTrace(lines, name, color, width, group, show=True, dash=None, hover=None):
        xs, ys, zs = [], [], []
        for ln in lines:
            xs += list(ln[:, 0]) + [None]
            ys += list(ln[:, 1]) + [None]
            zs += list(ln[:, 2]) + [None]
        return go.Scatter3d(x=xs, y=ys, z=zs, mode='lines', name=name, legendgroup=group,
                            line=dict(color=color, width=width, dash=dash), hoverinfo='skip' if hover is None else hover,
                            visible=True if show else 'legendonly')
    traces.append(polyTrace(cells, 'cell edges', '#c8c8c8', 1, 'grid'))
    traces.append(polyTrace(panels, 'panel edges', '#707070', 3, 'grid'))
    for arm in arms:
        traces.append(polyTrace([paths[arm][:, k] for k in range(len(x0))], ARM_LABELS[arm], COLORS[arm],
                                5 if arm == 'E' else 3, arm, dash='dash' if arm == 'E' else None))
        traces.append(go.Scatter3d(x=paths[arm][-1, :, 0], y=paths[arm][-1, :, 1], z=paths[arm][-1, :, 2],
                                   mode='markers', marker=dict(size=3, color=COLORS[arm], symbol='square'),
                                   name=f'{arm} end', legendgroup=arm, showlegend=False))
    traces.append(go.Scatter3d(x=x0[:NWORST, 0], y=x0[:NWORST, 1], z=x0[:NWORST, 2], mode='markers',
                               marker=dict(size=6, color='red', symbol='diamond'), name='start, worst Whitney paths'))
    traces.append(go.Scatter3d(x=x0[NWORST:, 0], y=x0[NWORST:, 1], z=x0[NWORST:, 2], mode='markers',
                               marker=dict(size=4, color='black'), name='start, random seeds'))
    fig3 = go.Figure(traces)
    fig3.update_layout(title=f'{args.case} trajectories, M={M} (click legend entries to toggle arms)',
                       scene=dict(aspectmode='data', xaxis_visible=False, yaxis_visible=False,
                                  zaxis_visible=False),
                       margin=dict(l=0, r=0, t=40, b=0), legend=dict(itemsizing='constant'))
    out = HERE / f'trajectories_{tag}.html'
    fig3.write_html(out, include_plotlyjs=True)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
