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

The module also has 4 more closed-loop tests (test_closed_loop_*) built
from the SAME loop shape, shifted to 4 different locations, to check that
the "closed loop -> 0" discrete-Stokes property survives the seams where
a bug would most likely hide: entirely inside one panel, straddling 2
panels, straddling 3 panels (at a cube vertex), and crossing the
longitude +-180 dateline. See the comment above LOOP_ONE_PANEL etc. for
the geometry.
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


# --- Closed loops that specifically exercise panel-boundary and dateline
# handling --------------------------------------------------------------
#
# The discrete-Stokes property (see module docstring) says a closed loop's
# line integral of the exact discrete gradient is ~0 regardless of
# resolution and regardless of which cell(s) it runs through. The 4 tests
# built from the loops below check that this actually holds not just for
# a loop entirely inside one cubed-sphere panel, but for the SAME loop
# shape shifted to straddle 2 panels, 3 panels, and one crossing the
# longitude +-180 dateline -- the seams where a bug in per-cell
# dateline/corner handling (see generate_cubedsphere_grid._fixCellLongitudes)
# would most likely show up.

def _faceOf(lon, lat):
    """Which of the 6 cube faces (+X, -X, +Y, -Y, +Z, -Z) a (lon, lat)
    direction falls on: whichever cube axis has the largest-magnitude
    component, with its sign -- the same classification implicit in
    generate_cubedsphere_grid._faceFrames."""
    lam, th = lon * DEG2RAD, lat * DEG2RAD
    d = numpy.array([numpy.cos(th) * numpy.cos(lam),
                     numpy.cos(th) * numpy.sin(lam),
                     numpy.sin(th)])
    axis = int(numpy.argmax(numpy.abs(d)))
    sign = '+' if d[axis] > 0. else '-'
    return sign + 'XYZ'[axis]


def _facesTouched(path, n=50):
    """The set of cube faces touched anywhere along a closed polyline's
    perimeter (not just its corners -- n interior samples per segment
    too), used to confirm each loop below actually straddles the panels
    it's meant to and no others."""
    faces = set()
    for (lon0, lat0, _), (lon1, lat1, _) in zip(path[:-1], path[1:]):
        for t in numpy.linspace(0., 1., n):
            faces.add(_faceOf(lon0 + t * (lon1 - lon0), lat0 + t * (lat1 - lat0)))
    return faces


def _closedRectangle(lon0, lat0, dlon, dlat):
    """A closed lon/lat rectangle centred at (lon0, lat0)."""
    return numpy.array([(lon0 - dlon, lat0 - dlat, 0.),
                        (lon0 + dlon, lat0 - dlat, 0.),
                        (lon0 + dlon, lat0 + dlat, 0.),
                        (lon0 - dlon, lat0 + dlat, 0.),
                        (lon0 - dlon, lat0 - dlat, 0.)])


_LOOP_HALFWIDTH = 8.0  # degrees, shared by the 3 panel-crossing loops below

# +X face is centred at (lon, lat) = (0, 0); this loop sits comfortably
# inside it, away from every panel boundary.
LOOP_ONE_PANEL = _closedRectangle(0., 0., _LOOP_HALFWIDTH, _LOOP_HALFWIDTH)

# +X and +Y share the meridian lon = 45 deg (for |lat| < 35.264 deg, see
# _CORNER_LAT below) -- the SAME loop, shifted to be centred on that
# meridian, straddles exactly those two panels.
LOOP_TWO_PANELS = _closedRectangle(45., 0., _LOOP_HALFWIDTH, _LOOP_HALFWIDTH)

# +X, +Y and +Z all meet at the cube vertex in direction (1, 1, 1)/sqrt(3),
# i.e. (lon, lat) = (45, 35.264 deg). A cube vertex always has exactly 3
# faces meeting there, so the SAME loop, shifted to be centred on this
# vertex, touches exactly those 3 panels and no 4th one.
_CORNER_LAT = numpy.degrees(numpy.arcsin(1. / numpy.sqrt(3.)))
LOOP_THREE_PANELS = _closedRectangle(45., _CORNER_LAT, _LOOP_HALFWIDTH, _LOOP_HALFWIDTH)

# -X face is centred at (lon, lat) = (180, 0) -- i.e. it's the panel the
# longitude +-180 branch cut runs straight through. This loop is fully
# inside -X (it stays well clear of the +Y/-Y panels, which start around
# lon = +-135), but its two vertical edges each cross the dateline once
# (e.g. from lon=170 to lon=190, the latter representing the same point as
# lon=-170).
LOOP_DATELINE = _closedRectangle(180., 0., 10., 10.)

# Confirm each loop touches exactly the panels it's meant to, at
# collection time -- if this ever fails, the geometry assumptions above
# (not the discrete-Stokes property under test) are what broke.
assert _facesTouched(LOOP_ONE_PANEL) == {'+X'}
assert _facesTouched(LOOP_TWO_PANELS) == {'+X', '+Y'}
assert _facesTouched(LOOP_THREE_PANELS) == {'+X', '+Y', '+Z'}
assert _facesTouched(LOOP_DATELINE) == {'-X'}


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


@pytest.mark.parametrize("M", [4, 8, 16])
def test_closed_loop_single_panel(M):
    """
    Test 1: closed loop entirely inside one panel (+X). Baseline
    discrete-Stokes check, with no panel boundary or dateline involved.
    """
    flux = _computeLineIntegral(M, LOOP_ONE_PANEL)
    print(f'M={M} single-panel loop: flux={flux:.3e} (exact 0)')
    assert abs(flux) < 1.e-8


@pytest.mark.parametrize("M", [4, 8, 16])
def test_closed_loop_two_panels(M):
    """
    Test 2: the same loop as test 1, shifted to straddle the +X/+Y panel
    boundary. Still exactly 0, i.e. the discrete gradient built
    independently on each panel is consistent across the shared edge.
    """
    flux = _computeLineIntegral(M, LOOP_TWO_PANELS)
    print(f'M={M} two-panel loop: flux={flux:.3e} (exact 0)')
    assert abs(flux) < 1.e-8


@pytest.mark.parametrize("M", [4, 8, 16])
def test_closed_loop_three_panels(M):
    """
    Test 3: the same loop again, shifted to the +X/+Y/+Z cube vertex, so
    it straddles all 3 panels that meet there.
    """
    flux = _computeLineIntegral(M, LOOP_THREE_PANELS)
    print(f'M={M} three-panel loop: flux={flux:.3e} (exact 0)')
    assert abs(flux) < 1.e-8


@pytest.mark.parametrize("M", [4, 8, 16])
def test_closed_loop_across_dateline(M):
    """
    Test 4: the same loop shape, now centred on (lon, lat) = (180, 0) --
    fully inside the -X panel, but with both its vertical edges crossing
    the longitude +-180 branch cut once each.
    """
    flux = _computeLineIntegral(M, LOOP_DATELINE)
    print(f'M={M} dateline loop: flux={flux:.3e} (exact 0)')
    assert abs(flux) < 1.e-8


if __name__ == '__main__':
    test_open_path_crossing_lon0()
    for M in (8, 16, 32, 64):
        test_open_path_exact_when_endpoints_are_grid_nodes(M)
    for M in (4, 8, 16):
        test_closed_path_crossing_lon0(M)
        test_closed_loop_single_panel(M)
        test_closed_loop_two_panels(M)
        test_closed_loop_three_panels(M)
        test_closed_loop_across_dateline(M)
