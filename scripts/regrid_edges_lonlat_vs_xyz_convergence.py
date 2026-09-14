#!/usr/bin/env python
"""
Compares mint.RegridEdges convergence -- both ACCURACY and WALL-CLOCK TIME --
between the lon-lat locator (useXYZLocator=False, the original code path) and
the genuine xyz locator (useXYZLocator=True, see
project_mint_locator_3d_fix.md) on the SAME cubed-sphere panel, for the SAME
physical field, expressed consistently in each representation's own raw
coordinate basis.

Two panels are compared:
- '+Z': pole-containing. mint's earlier VectorInterp/PolylineIntegral
  convergence studies (project_mint_whitney_convergence_studies.md) found
  lon-lat reconstruction degrades to ~O(h^1.9) here (a genuine coordinate-
  basis-distortion artifact near the pole, not a modelling limitation) --
  this script checks whether the SAME degradation shows up in RegridEdges
  (not just pointwise VectorInterp reconstruction), and whether the xyz
  locator avoids it.
- '+X': equatorial, no coordinate singularity in either representation --
  expected control case where lon-lat and xyz should track each other
  closely, isolating the pole panel's result as the real effect rather than
  some other lon-lat-vs-xyz artifact.

For each panel, at each destination resolution Mdst, a source grid 4x finer
(Msrc = 4*Mdst) is regridded onto the destination grid -- the same "coarse
destination, fine source" scenario scripts/xyz_regrid_prototype.py validated
against an exact analytic reference, now using the real, compiled
mint.RegridEdges (both locator paths) rather than a from-scratch Python
prototype. RMS/max edge-circulation error against the exact field is the
accuracy metric; wall-clock time for buildLocator and computeWeights
(median of a few repeats) is the timing metric -- the natural follow-up
question once both paths are real, tested code rather than one prototype
and one library.

Field: same pole-asymmetric vector potential u = curl(psi * r_hat),
psi(lambda, theta) = cos(theta)*(1+sin(theta))*cos(lambda), used throughout
this project's convergence studies (pole_asymmetric_convergence.py); reused
here via vector_potential_reconstruction_convergence.exactVector (lon-lat,
DEG2RAD-corrected raw-degree-basis components) and
xyz_pole_panel_validation.vectorFieldXYZ (genuine ambient xyz components) --
the SAME physical field in the two representations, so the comparison is
apples-to-apples.

Usage: python scripts/regrid_edges_lonlat_vs_xyz_convergence.py
"""
import sys
import time
from pathlib import Path

import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mint import Grid, RegridEdges, CELL_BY_CELL_DATA

HERE = Path(__file__).absolute().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from whitney2 import edgeIntegral  # noqa: E402
from generate_cubedsphere_grid import _faceFrames, _fixCellLongitudes  # noqa: E402
from vector_potential_reconstruction_convergence import exactVector  # noqa: E402
from xyz_pole_panel_validation import vectorFieldXYZ  # noqa: E402

NUM_VERTS_PER_QUAD = 4
# same (normal, u_hat, v_hat) face frames generate_cubedsphere_grid uses for
# the whole-sphere lon-lat build -- reused directly (not re-derived) so the
# lon-lat and xyz panel builders below are GUARANTEED to describe the exact
# same physical patch of the sphere, corner for corner.
_FRAMES0 = {name: (normal, u_hat, v_hat) for name, normal, u_hat, v_hat in _faceFrames(0.0)}


def vectorFieldLonLat(pts):
    """(..., 2) (lon_deg, lat_deg) -> (..., 2) raw-coordinate-basis
    (u_lam, u_theta) -- same physical field as vectorFieldXYZ, reusing the
    already DEG2RAD-corrected exactVector."""
    v3 = exactVector(pts[..., 0], pts[..., 1])
    return v3[..., :2]


def _projectPanel(M, panel):
    """Shared geometry step for both builders below: M+1 x M+1 gnomonic
    projection of one cube face onto the unit sphere."""
    normal, u_hat, v_hat = _FRAMES0[panel]
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = numpy.meshgrid(edges, edges, indexing='ij')
    cube_pts = (normal[None, None, :] + U[..., None] * u_hat[None, None, :]
                + V[..., None] * v_hat[None, None, :])
    norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
    return cube_pts / norm  # (M+1, M+1, 3) unit-sphere xyz


def buildPanelPointsLonLat(M, panel):
    """(M*M, 4, 3) corners: (lon_deg, lat_deg, 0), with the same per-cell
    dateline/pole fix-up generate_cubedsphere_grid.cubedSphereGridPoints
    applies (reused directly, not reimplemented)."""
    sphere_pts = _projectPanel(M, panel)
    lon = numpy.degrees(numpy.arctan2(sphere_pts[..., 1], sphere_pts[..., 0]))
    lat = numpy.degrees(numpy.arcsin(numpy.clip(sphere_pts[..., 2], -1.0, 1.0)))
    points = numpy.empty((M * M, NUM_VERTS_PER_QUAD, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            corner_lon = numpy.array([lon[i, j], lon[i + 1, j], lon[i + 1, j + 1], lon[i, j + 1]])
            corner_lat = numpy.array([lat[i, j], lat[i + 1, j], lat[i + 1, j + 1], lat[i, j + 1]])
            points[icell, :, 0] = _fixCellLongitudes(corner_lon, corner_lat)
            points[icell, :, 1] = corner_lat
            points[icell, :, 2] = 0.0
            icell += 1
    return points


def buildPanelPointsXYZ(M, panel):
    """(M*M, 4, 3) corners: genuine unit-sphere (x, y, z), no lon/lat at
    all -- same panel/resolution as buildPanelPointsLonLat above."""
    sphere_pts = _projectPanel(M, panel)
    points = numpy.empty((M * M, NUM_VERTS_PER_QUAD, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            points[icell, 0] = sphere_pts[i, j]
            points[icell, 1] = sphere_pts[i + 1, j]
            points[icell, 2] = sphere_pts[i + 1, j + 1]
            points[icell, 3] = sphere_pts[i, j + 1]
            icell += 1
    return points


def buildGridLonLat(points):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1, degrees=True)
    grid.setPoints(points)
    return grid


def buildGridXYZ(points):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=0)
    grid.setPoints(points)
    return grid


def buildExactEdgeData(points, vectorField, ncomp, n=16):
    """Cell-by-cell circulations, straight-chord-in-this-representation
    convention (matches what RegridEdges itself treats as "the edge") --
    edges 2 and 3 canonicalized to the positive-parametric-direction sign
    RegridEdges' own mntQuadEdgeIter.cpp uses (3-->2, 0-->3), same
    convention as mint/tests/test_regrid_edges_xyz.py's buildExactEdgeData."""
    pts = points[..., :ncomp]
    ncells = pts.shape[0]
    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        for icell in range(ncells):
            data[icell, i0] = sign * edgeIntegral(pts[icell, i0], pts[icell, i1], vectorField, None, n)
    return data


# (label, buildPoints, buildGrid, vectorField, ncomp)
COORD_SYSTEMS = {
    'lon-lat': (buildPanelPointsLonLat, buildGridLonLat, vectorFieldLonLat, 2),
    'xyz':     (buildPanelPointsXYZ, buildGridXYZ, vectorFieldXYZ, 3),
}


def regridOnce(coordName, panel, Msrc, Mdst, numRepeats=3):
    """One src->dst regrid in one coordinate representation: returns
    (rms, max, tBuildLocator, tComputeWeights, tApply) -- the three timings
    are the MEDIAN of numRepeats independent buildLocator+computeWeights+
    apply runs (a fresh RegridEdges/grid each time, to not favour either
    path with warm caches from a previous repeat)."""
    buildPoints, buildGrid, vectorField, ncomp = COORD_SYSTEMS[coordName]
    useXYZLocator = (coordName == 'xyz')

    srcPoints = buildPoints(Msrc, panel)
    dstPoints = buildPoints(Mdst, panel)
    srcData = buildExactEdgeData(srcPoints, vectorField, ncomp).flatten()
    exact = buildExactEdgeData(dstPoints, vectorField, ncomp)

    tBuild, tWeights, tApply = [], [], []
    rms = mx = None
    for _ in range(numRepeats):
        srcGrid = buildGrid(srcPoints)
        dstGrid = buildGrid(dstPoints)

        rg = RegridEdges()
        rg.setSrcGrid(srcGrid)
        rg.setDstGrid(dstGrid)

        t0 = time.perf_counter()
        rg.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=0,
                         useXYZLocator=useXYZLocator)
        tBuild.append(time.perf_counter() - t0)

        t0 = time.perf_counter()
        rg.computeWeights()
        tWeights.append(time.perf_counter() - t0)

        dstData = numpy.zeros(Mdst * Mdst * NUM_VERTS_PER_QUAD)
        t0 = time.perf_counter()
        rg.apply(srcData, dstData, placement=CELL_BY_CELL_DATA)
        tApply.append(time.perf_counter() - t0)

        dstData = dstData.reshape(Mdst * Mdst, NUM_VERTS_PER_QUAD)
        err = dstData - exact
        rms = float(numpy.sqrt(numpy.mean(err ** 2)))
        mx = float(numpy.max(numpy.abs(err)))

    return rms, mx, float(numpy.median(tBuild)), float(numpy.median(tWeights)), float(numpy.median(tApply))


def fitAlpha(Ms, errs):
    """Least-squares fit of log(err) = -alpha*log(M) + c."""
    Ms = numpy.asarray(Ms, dtype=numpy.float64)
    errs = numpy.asarray(errs, dtype=numpy.float64)
    slope, intercept = numpy.polyfit(numpy.log(Ms), numpy.log(errs), 1)
    return -slope, numpy.exp(intercept)


PANELS = ('+Z', '+X')
RESOLUTIONS_MDST = (2, 4, 8, 16)  # Msrc = 4 * Mdst for each


def main():
    results = {}  # (panel, coordName) -> dict of lists
    for panel in PANELS:
        for coordName in COORD_SYSTEMS:
            r = {'Mdst': [], 'Msrc': [], 'rms': [], 'max': [],
                 'tBuild': [], 'tWeights': [], 'tApply': []}
            for Mdst in RESOLUTIONS_MDST:
                Msrc = 4 * Mdst
                rms, mx, tBuild, tWeights, tApply = regridOnce(coordName, panel, Msrc, Mdst)
                r['Mdst'].append(Mdst)
                r['Msrc'].append(Msrc)
                r['rms'].append(rms)
                r['max'].append(mx)
                r['tBuild'].append(tBuild)
                r['tWeights'].append(tWeights)
                r['tApply'].append(tApply)
                print(f'panel={panel:2s} coords={coordName:7s} Msrc={Msrc:3d} Mdst={Mdst:3d}  '
                      f'rms={rms:.3e}  max={mx:.3e}  '
                      f'tBuild={tBuild*1e3:7.2f}ms  tWeights={tWeights*1e3:7.2f}ms  tApply={tApply*1e3:6.2f}ms')
            results[(panel, coordName)] = r

    print()
    for panel in PANELS:
        for coordName in COORD_SYSTEMS:
            r = results[(panel, coordName)]
            alpha, _ = fitAlpha(r['Mdst'], r['rms'])
            totalTime = numpy.array(r['tBuild']) + numpy.array(r['tWeights'])
            print(f'panel={panel:2s} coords={coordName:7s}: RMS convergence order alpha={alpha:.2f}  '
                  f'(build+computeWeights time at largest Mdst={r["Mdst"][-1]}: {totalTime[-1]*1e3:.1f}ms)')

    # --- plots: 2 rows (accuracy, timing) x 2 columns (one per panel) -----
    fig, axes = plt.subplots(2, len(PANELS), figsize=(6 * len(PANELS), 9))
    colors = {'lon-lat': 'C0', 'xyz': 'C1'}
    for col, panel in enumerate(PANELS):
        axAcc = axes[0, col]
        axTime = axes[1, col]
        for coordName in COORD_SYSTEMS:
            r = results[(panel, coordName)]
            alpha, _ = fitAlpha(r['Mdst'], r['rms'])
            axAcc.loglog(r['Mdst'], r['rms'], 'o-', color=colors[coordName],
                         label=f'{coordName} (alpha={alpha:.2f})')
            totalTime = numpy.array(r['tBuild']) + numpy.array(r['tWeights'])
            axTime.loglog(r['Mdst'], totalTime * 1e3, 's-', color=colors[coordName], label=coordName)
        # reference slopes on the accuracy plot
        Ms = numpy.array(RESOLUTIONS_MDST, dtype=numpy.float64)
        ref = results[(panel, 'lon-lat')]['rms'][0] * (Ms[0] / Ms) ** 2
        axAcc.loglog(Ms, ref, 'k--', alpha=0.5, label='$M^{-2}$ reference')

        axAcc.set_title(f"panel {panel}{'  (pole-containing)' if panel == '+Z' else '  (equatorial, control)'}")
        axAcc.set_xlabel('Mdst (Msrc = 4*Mdst)')
        axAcc.set_ylabel('RMS edge-circulation error')
        axAcc.legend()
        axAcc.grid(True, which='both', alpha=0.3)

        axTime.set_xlabel('Mdst (Msrc = 4*Mdst)')
        axTime.set_ylabel('buildLocator + computeWeights [ms]')
        axTime.legend()
        axTime.grid(True, which='both', alpha=0.3)

    fig.suptitle('RegridEdges: lon-lat vs xyz locator, accuracy and timing')
    fig.tight_layout()
    outPng = HERE / 'regrid_edges_lonlat_vs_xyz_convergence.png'
    fig.savefig(outPng, dpi=120)
    print(f'\nwrote {outPng}')


if __name__ == '__main__':
    main()
