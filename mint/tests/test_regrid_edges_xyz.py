"""
End-to-end regression test for RegridEdges with useXYZLocator=True -- the
real C++ pipeline (RegridEdges -> PolysegmentIter -> vmtXYZCellLocator::
findIntersectionsWithLine), not the Python prototype
(scripts/xyz_regrid_prototype.py) it was validated against first.

Regrids a genuinely (x, y, z) embedded cubed-sphere '+X' panel from a finer
source resolution to a coarser destination resolution, and checks the
regridded destination edge circulations against the exact analytic values.

Tolerance: NOT machine precision. Per the design discussion (and the
prototype's own validation), treating a destination edge as a straight 3D
chord introduces a real, expected, resolution-dependent approximation error
(the chord dips slightly inside the curved surface except at its own
endpoints) -- the prototype measured RMS ~1e-2 to ~1e-3 at comparable
resolutions, decreasing at roughly O(h^2.5-3.5) as destination resolution
increases, and only breaking down at deliberately unrealistic resolutions
(e.g. one destination cell spanning most of a 90-degree panel). This test
picks a tolerance well inside that envelope, not zero.
"""
import sys
from pathlib import Path

import numpy

from mint import Grid, RegridEdges, CELL_BY_CELL_DATA

HERE = Path(__file__).absolute().parent
SCRIPTS = HERE.parent.parent / 'scripts'
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
from whitney2 import edgeIntegral  # noqa: E402
from xyz_pole_panel_validation import vectorFieldXYZ  # noqa: E402


def buildPlusXPanelXYZ(M):
    """M x M cells covering the '+X' cube face, gnomonically projected onto
    the unit sphere -- same construction used throughout this codebase's
    xyz validation scripts/tests."""
    normal, u_hat, v_hat = numpy.array([1., 0., 0.]), numpy.array([0., 1., 0.]), numpy.array([0., 0., 1.])
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = numpy.meshgrid(edges, edges, indexing='ij')
    cube_pts = (normal[None, None, :] + U[..., None] * u_hat[None, None, :]
                + V[..., None] * v_hat[None, None, :])
    norm = numpy.linalg.norm(cube_pts, axis=-1, keepdims=True)
    sphere_pts = cube_pts / norm

    points = numpy.empty((M * M, 4, 3))
    icell = 0
    for i in range(M):
        for j in range(M):
            points[icell, 0] = sphere_pts[i, j]
            points[icell, 1] = sphere_pts[i + 1, j]
            points[icell, 2] = sphere_pts[i + 1, j + 1]
            points[icell, 3] = sphere_pts[i, j + 1]
            icell += 1
    return points


def buildXYZGrid(points):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=0)
    grid.setPoints(points)
    return grid


def buildExactEdgeData(points, n=16):
    """Cell-by-cell circulations in RegridEdges' own edge convention: edges
    are always canonicalized to point in the *positive parametric
    direction* (mntQuadEdgeIter.cpp), i.e. 0-->1, 1-->2, 3-->2, 0-->3 --
    NOT the raw corner_i -> corner_{i+1} cycle (which would have edges 2
    and 3 backwards, 2-->3 and 3-->0). Same sign-flip convention as
    test_vector_interp_xyz.py's constantFieldEdgeData."""
    ncells = points.shape[0]
    data = numpy.zeros((ncells, 4))
    for i0 in range(4):
        i1 = (i0 + 1) % 4
        sign = 1 - 2 * (i0 // 2)  # +1 for i0=0,1; -1 for i0=2,3
        for icell in range(ncells):
            data[icell, i0] = sign * edgeIntegral(points[icell, i0], points[icell, i1], vectorFieldXYZ, None, n)
    return data


def test_regrid_edges_xyz_cubedsphere_panel():
    Msrc, Mdst = 8, 4

    srcPoints = buildPlusXPanelXYZ(Msrc)
    dstPoints = buildPlusXPanelXYZ(Mdst)

    srcGrid = buildXYZGrid(srcPoints)
    dstGrid = buildXYZGrid(dstPoints)

    rg = RegridEdges()
    rg.setSrcGrid(srcGrid)
    rg.setDstGrid(dstGrid)
    rg.buildLocator(numCellsPerBucket=128, useXYZLocator=True)
    rg.computeWeights()

    srcData = buildExactEdgeData(srcPoints).flatten()
    dstData = numpy.zeros(Mdst * Mdst * 4)
    rg.apply(srcData, dstData, placement=CELL_BY_CELL_DATA)
    dstData = dstData.reshape(Mdst * Mdst, 4)

    exact = buildExactEdgeData(dstPoints)

    err = dstData - exact
    rms = float(numpy.sqrt(numpy.mean(err ** 2)))
    mx = float(numpy.max(numpy.abs(err)))

    # generous relative to the field's own magnitude (~O(1) circulations
    # per edge at this resolution) -- see module docstring for why this
    # isn't machine precision
    assert rms < 0.05, f'RMS regridding error too large: {rms:.3e}'
    assert mx < 0.15, f'max regridding error too large: {mx:.3e}'


if __name__ == '__main__':
    test_regrid_edges_xyz_cubedsphere_panel()
