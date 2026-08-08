"""
Line integral of an exact potential across the longitude = 0 line, on the
top (+Z) panel of a cubed-sphere grid.

The potential is

    V(lon, lat) = sin(lon) * cos(lat)**2

and the edge data fed to mint.PolylineIntegral is built from the EXACT
finite differences of V at each cell's 4 corners (see _computeEdgeData
below). mint's polyline integral is a mimetic (Whitney-edge-element)
construction, which has one exactness property regardless of resolution
or of how nonlinear the underlying potential is: the discrete curl of a
discrete gradient is exactly zero (discrete Stokes theorem). That gives
two different, and differently strict, expectations for the two paths
tested here:

  * CLOSED path crossing longitude = 0 (twice, net zero winding, and not
    encircling the pole): integral == 0 to near machine precision, at
    ANY resolution -- this is the discrete-Stokes property directly.
  * OPEN path crossing longitude = 0 once, with ARBITRARY endpoints (not
    grid vertices): integral -> V(last point) - V(first point) as the
    grid is refined, but only approximately (2nd order in 1/M) at finite
    resolution. Reconstructing the field *inside* a cell from its 4
    corner values is exact along that cell's own edges (an edge's data
    literally IS the exact potential difference between its two
    endpoints, for ANY potential, by the fundamental theorem of
    calculus), but an arbitrary interior point needs an interpolated
    value, which is only exact for a field that is linear in the cell's
    local coordinates -- sin(lon)*cos(lat)^2 is not. That leftover
    interpolation error sits at the path's two endpoints; it doesn't
    cancel because the two endpoints are different points with (in
    general) different leftover errors.
  * OPEN path crossing longitude = 0 once, with BOTH endpoints exactly
    ON grid vertices: integral == V(last point) - V(first point) to
    near machine precision, at ANY resolution, REGARDLESS of how the
    interior of the path cuts through cells (it need not follow mesh
    edges at all). No interpolation is needed at either endpoint --
    each is a grid node, and the per-edge data are already exact node
    potential differences there -- so the only source of error above
    disappears identically. (Contrast with mint/tests/test_polyline_integral.py's
    test_simple, which uses a genuinely linear potential and so gets an
    exact match for an open path even with non-vertex endpoints.)

Both paths are chosen to sit well inside the top panel. The top face
(normal +Z, see generate_cubedsphere_grid._faceFrames) covers directions
with tan(lat) >= max(|cos(lon)|, |sin(lon)|); at lon = 0 that boundary is
lat = 45 deg, tightening to lat = 35.264 deg at the face's corners
(lon = +-45, +-135). Both paths below stay at |lon| <= 30 and
lat >= 55 deg, comfortably clear of that boundary (worst case, at
lon = +-30, the boundary is at lat = 40.9 deg).
"""
import sys
from pathlib import Path

import numpy
import pytest

from mint import PolylineIntegral, CELL_BY_CELL_DATA

# scripts/generate_cubedsphere_grid.py is not part of the mint package;
# reach it via sys.path, the same way scripts/pole_asymmetric_paths_3d.py
# reaches its sibling scripts.
SCRIPTS_DIR = Path(__file__).absolute().parents[2] / 'scripts'
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
from generate_cubedsphere_grid import buildCubedSphereGrid, NUM_VERTS_PER_QUAD  # noqa: E402

DEG2RAD = numpy.pi / 180.0

# Open path: -30 -> 0 -> +30 deg longitude at high latitude, apex poleward
# of its endpoints. Crosses longitude = 0 once.
OPEN_PATH = numpy.array([(-30., 55., 0.),
                         (0., 70., 0.),
                         (30., 55., 0.)])

# Closed path: a lon/lat rectangle straddling longitude = 0. Its two
# horizontal edges each cross longitude = 0 (in opposite directions), for
# a net winding of zero and no pole enclosed.
CLOSED_PATH = numpy.array([(-30., 55., 0.),
                           (30., 55., 0.),
                           (30., 75., 0.),
                           (-30., 75., 0.),
                           (-30., 55., 0.)])

# For the "endpoints exactly on grid vertices" test: index pair (i, j) into
# the top face's own (M+1) x (M+1) node grid, expressed as fixed fractions
# of M so the same fractions locate a valid pair of nodes at any M
# divisible by 4. i is held fixed (an interior grid line, comfortably off
# the face boundary) while j is offset to either side of the row where
# V = 0 (longitude = 0), so the two nodes straddle longitude = 0.
_NODE_I_FRAC = 3. / 4.
_NODE_J_NEG_FRAC = 1. / 4.
_NODE_J_POS_FRAC = 3. / 4.


def _topFaceNode(M, i, j):
    """
    The (lon, lat) of node (i, j) of the top (+Z) face's own (M+1) x (M+1)
    grid -- i.e. the SAME formula generate_cubedsphere_grid.cubedSphereGridPoints
    uses for the '+Z' face (normal Z, u_hat=X, v_hat=Y), so this returns a
    point that is bit-identical to an actual vertex of the mesh built by
    buildCubedSphereGrid(M). That bit-identity, not just numerical
    closeness, is what makes the exactness test below meaningful.
    """
    edges = numpy.linspace(-1., 1., M + 1)
    U, V = edges[i], edges[j]
    pt = numpy.array([U, V, 1.0])
    pt = pt / numpy.linalg.norm(pt)
    lon = numpy.degrees(numpy.arctan2(pt[1], pt[0]))
    lat = numpy.degrees(numpy.arcsin(pt[2]))
    return lon, lat


def _onNodeOpenPath(M):
    """Open, node-to-node path straddling longitude = 0, cutting diagonally
    through whatever cells lie between the two nodes (deliberately NOT
    following mesh edges)."""
    i = round(_NODE_I_FRAC * M)
    j_neg = round(_NODE_J_NEG_FRAC * M)
    j_pos = round(_NODE_J_POS_FRAC * M)
    lonA, latA = _topFaceNode(M, i, j_neg)
    lonB, latB = _topFaceNode(M, i, j_pos)
    assert lonA < 0. < lonB, 'expected the two nodes to straddle longitude = 0'
    return numpy.array([(lonA, latA, 0.), (lonB, latB, 0.)])


def potential(p):
    """Exact global potential V(lon, lat) = sin(lon) * cos(lat)**2."""
    lam = p[..., 0] * DEG2RAD
    theta = p[..., 1] * DEG2RAD
    return numpy.sin(lam) * numpy.cos(theta) ** 2


def _computeEdgeData(points):
    """
    Build CELL_BY_CELL_DATA edge data from the exact potential, following
    the node/edge ordering documented in test_polyline_integral.py:

        3-->--2
        |     |
        ^     ^
        |     |
        0-->--1

    i.e. edges 0 and 1 run "forward" (node i -> node i+1) while edges 2
    and 3 run "backward" (node i+1 -> node i), hence the sign flip for
    i0 >= 2.
    """
    ncells = points.shape[0]
    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * (potential(points[:, i1, :]) - potential(points[:, i0, :]))
    return data


def _computeLineIntegral(M, path):
    grid = buildCubedSphereGrid(M)
    points = grid.getPoints()
    data = _computeEdgeData(points)

    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=False)
    pli.computeWeights(path, counterclock=False)
    return pli.getIntegral(data, placement=CELL_BY_CELL_DATA)


def test_open_path_crossing_lon0():
    """
    Open path from (-30, 55) to (30, 55), crossing longitude = 0 at its
    apex (0, 70). The computed line integral should converge to the
    exact V(end) - V(start) as the grid is refined (2nd order in 1/M,
    since the potential is not linear in each cell's local coordinates
    -- see module docstring), and the error should shrink monotonically
    with M.
    """
    exact = potential(OPEN_PATH[-1]) - potential(OPEN_PATH[0])
    assert abs(exact) > 1.e-3, 'expected a nonzero exact flux for this test to be meaningful'

    errors = []
    for M in (4, 8, 16, 32):
        flux = _computeLineIntegral(M, OPEN_PATH)
        err = abs(flux - exact)
        tol = 2.0 / M ** 2
        print(f'M={M} open path: flux={flux:.9f} exact={exact:.9f} '
              f'error={err:.3e} (tol {tol:.3e})')
        assert err < tol, f'error {err:.3e} exceeds the O(1/M^2) convergence tolerance {tol:.3e}'
        errors.append(err)

    for M_prev, M_next, e_prev, e_next in zip((4, 8, 16), (8, 16, 32), errors[:-1], errors[1:]):
        assert e_next < e_prev, (
            f'refining M={M_prev} -> M={M_next} should reduce the error '
            f'({e_prev:.3e} -> {e_next:.3e}): the scheme is not converging')


@pytest.mark.parametrize("M", [8, 16, 32, 64])
def test_open_path_exact_when_endpoints_are_grid_nodes(M):
    """
    Same kind of open path as test_open_path_crossing_lon0 (crosses
    longitude = 0 once), but with both endpoints placed exactly on grid
    vertices instead of arbitrary points. Unlike that test, the integral
    should now match the exact V(end) - V(start) to near machine
    precision at EVERY resolution, not just in the limit -- see the
    module docstring for why.
    """
    path = _onNodeOpenPath(M)
    exact = potential(path[-1]) - potential(path[0])
    assert abs(exact) > 1.e-3, 'expected a nonzero exact flux for this test to be meaningful'

    flux = _computeLineIntegral(M, path)
    err = abs(flux - exact)
    print(f'M={M} on-node open path: flux={flux:.12f} exact={exact:.12f} error={err:.3e}')
    assert err < 1.e-8, f'error {err:.3e} is not near machine precision: endpoints may not be exact grid nodes'


@pytest.mark.parametrize("M", [4, 8, 16])
def test_closed_path_crossing_lon0(M):
    """
    Closed rectangular path straddling longitude = 0 (crossed twice, net
    zero winding, pole not enclosed). The line integral of an exact field
    around any closed loop is zero.
    """
    flux = _computeLineIntegral(M, CLOSED_PATH)
    print(f'M={M} closed path: flux={flux:.3e} (exact 0)')
    assert abs(flux) < 1.e-8


if __name__ == '__main__':
    test_open_path_crossing_lon0()
    for M in (8, 16, 32, 64):
        test_open_path_exact_when_endpoints_are_grid_nodes(M)
    for M in (4, 8, 16):
        test_closed_path_crossing_lon0(M)
