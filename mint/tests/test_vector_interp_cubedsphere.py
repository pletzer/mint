"""
Reconstruct the vector field grad(chi) on a cubed-sphere grid, from its
edge-integral (1-form) representation, using mint.VectorInterp -- and
check the reconstruction converges as the grid is refined.

This continues the same 0-form/1-form/vector story as the other
cubed-sphere tests in this directory (see test_cubedsphere.py):

  * chi(lon, lat) = sin(lon) * cos(lat)**2 is the 0-form (potential),
    evaluated at grid nodes.
  * The edge data (see _computeEdgeData) is the exact 1-form dchi,
    i.e. the same "exact finite differences of chi at each cell's 4
    corners" used throughout test_cubedsphere.py: data[edge] is exactly
    chi(edge's end node) - chi(edge's start node), for ANY potential, by
    the fundamental theorem of calculus.
  * mint.VectorInterp.getEdgeVectors reconstructs a VECTOR from that
    1-form's edge data (Whitney-1-form basis functions, see the
    conversation with Claude on duality/Kronecker-delta degrees of
    freedom): the discrete analogue of "raise the index" on dchi to get
    grad(chi).

Three-part structure:

  1. exactGradient(lonDeg, latDeg) computes the EXACT vector field
     components: the partial derivatives of chi with respect to the
     grid's own (lon, lat) coordinates (in degrees, matching
     mint.Grid's `degrees=True` flag), i.e. grad(chi) expressed in the
     grid's native coordinate basis -- (dchi/dlon, dchi/dlat, 0). This
     is deliberately NOT converted to physical zonal/meridional wind
     (which would need an extra 1/cos(lat) metric correction, see
     mint/tests/test_vector_interp.py's test_accuracy): mint.VectorInterp
     itself only ever sees the grid's raw (lon, lat, 0) point
     coordinates and the edge data's raw numbers, with no notion of a
     spherical metric, so THIS is the fair, apples-to-apples target for
     what it reconstructs. It's also exactly what the edge data was
     built to be consistent with, so a path integral of this same exact
     vector field along any mesh edge reproduces that edge's data value
     precisely (as in test_cubedsphere.py).
  2. _reconstructVectors(M) builds the grid, the edge data, and calls
     mint.VectorInterp.getEdgeVectors to reconstruct the vector at each
     cell's own centre (a target point strictly inside that cell, at
     parametric (xi, eta) = (0.5, 0.5), i.e. the plain average of its 4
     corners).
  3. test_convergence_with_resolution checks the error between (1) and
     (2) shrinks as the cubed-sphere is refined (M = 4, 8, 16, 32).

Whitney-1-form reconstruction inside a cell is exact only for a field
that's bilinear in the cell's local coordinates (see
test_cubedsphere.py's test_open_path_crossing_lon0 docstring) --
sin(lon)*cos(lat)^2 is not, so the reconstruction error is only expected
to shrink with resolution, not vanish outright. It does so at two
different rates:

  * RMS error over all cells shrinks like 1/M^2 (2nd order): the typical
    cell is only mildly distorted, so ordinary bilinear-interpolation
    error theory applies.
  * MAX error only shrinks like 1/M (1st order), because the very worst
    cells are always the ones immediately surrounding the +-Z panels'
    pole vertex: (lon, lat) itself is a singular coordinate there
    (lines of constant longitude converge on the pole), which distorts
    the *coordinate representation* of grad(chi) far more than ordinary
    cell curvature does elsewhere on the sphere, independent of how
    fine the grid is. This was confirmed directly: at M=16 the 5 worst
    cells are all within ~6 deg of a pole.

Both rates are still genuine convergence (the error -> 0 as M -> inf),
just at different orders -- the tolerances below reflect that (1/M^2 for
the RMS check, 1/M for the max check), rather than using one common,
overly loose bound for both.
"""
import sys
from pathlib import Path

import numpy
import pytest

from mint import VectorInterp, CELL_BY_CELL_DATA

# scripts/generate_cubedsphere_grid.py is not part of the mint package;
# reach it via sys.path, as in test_cubedsphere.py.
SCRIPTS_DIR = Path(__file__).absolute().parents[2] / 'scripts'
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
from generate_cubedsphere_grid import buildCubedSphereGrid, NUM_VERTS_PER_QUAD  # noqa: E402

DEG2RAD = numpy.pi / 180.0


def potential(p):
    """The 0-form: chi(lon, lat) = sin(lon) * cos(lat)**2."""
    lam = p[..., 0] * DEG2RAD
    theta = p[..., 1] * DEG2RAD
    return numpy.sin(lam) * numpy.cos(theta) ** 2


def exactGradient(lonDeg, latDeg):
    """
    Part (1): the exact vector field, i.e. grad(chi) expressed in the
    grid's own (lon_deg, lat_deg) coordinate basis:

        dchi/dlon_deg =  DEG2RAD * cos(lon*DEG2RAD) * cos(lat*DEG2RAD)**2
        dchi/dlat_deg = -DEG2RAD * sin(lon*DEG2RAD) * sin(2*lat*DEG2RAD)

    (the DEG2RAD factors are the chain rule for differentiating with
    respect to the degree-valued coordinate rather than radians; the
    second uses sin(2x) = 2 sin(x) cos(x)). Accepts scalar or array
    lonDeg/latDeg; returns an array of shape (..., 3) with a trailing
    zero z-component, matching mint.VectorInterp's output shape.
    """
    lam = lonDeg * DEG2RAD
    theta = latDeg * DEG2RAD
    dlon = DEG2RAD * numpy.cos(lam) * numpy.cos(theta) ** 2
    dlat = -DEG2RAD * numpy.sin(lam) * numpy.sin(2. * theta)
    return numpy.stack([dlon, dlat, numpy.zeros_like(dlon)], axis=-1)


def test_exact_gradient_matches_finite_difference():
    """
    Sanity-check exactGradient itself (independent of the grid/VectorInterp
    machinery below) against a centred finite difference of potential, at
    an arbitrary, non-critical point.
    """
    lon0, lat0, eps = 37.0, 21.0, 1.e-4
    fd_dlon = (potential(numpy.array([lon0 + eps, lat0, 0.]))
               - potential(numpy.array([lon0 - eps, lat0, 0.]))) / (2. * eps)
    fd_dlat = (potential(numpy.array([lon0, lat0 + eps, 0.]))
               - potential(numpy.array([lon0, lat0 - eps, 0.]))) / (2. * eps)
    exact = exactGradient(lon0, lat0)
    print(f'finite-diff grad = ({fd_dlon:.9f}, {fd_dlat:.9f}), '
          f'exactGradient = ({exact[0]:.9f}, {exact[1]:.9f})')
    assert abs(exact[0] - fd_dlon) < 1.e-6
    assert abs(exact[1] - fd_dlat) < 1.e-6


def _computeEdgeData(points):
    """
    The 1-form dchi as CELL_BY_CELL_DATA, following the same node/edge
    ordering as test_cubedsphere.py's _computeEdgeData:

        3-->--2
        |     |
        ^     ^
        |     |
        0-->--1
    """
    ncells = points.shape[0]
    data = numpy.zeros((ncells, NUM_VERTS_PER_QUAD))
    for i0 in range(NUM_VERTS_PER_QUAD):
        i1 = (i0 + 1) % NUM_VERTS_PER_QUAD
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * (potential(points[:, i1, :]) - potential(points[:, i0, :]))
    return data


def _cellCenters(points):
    """Each cell's own centre: the plain average of its 4 corners, i.e.
    parametric (xi, eta) = (0.5, 0.5) -- always strictly inside the cell
    (never on a shared edge/vertex), and safe to average directly since
    each cell's own corner longitudes are already made self-consistent
    by generate_cubedsphere_grid._fixCellLongitudes."""
    return points[:, :, :2].mean(axis=1)  # (ncells, 2): (lon, lat)


def _reconstructVectors(M):
    """
    Part (2): build the grid + exact 1-form edge data, then use
    mint.VectorInterp to reconstruct the vector at each cell's own
    centre.

    :returns: (centers, numeric, exact) -- each an (ncells, ...) array:
              the target points' (lon, lat), VectorInterp's output, and
              exactGradient evaluated at the same points.
    """
    grid = buildCubedSphereGrid(M)
    points = grid.getPoints()
    data = _computeEdgeData(points)
    centers = _cellCenters(points)

    targets = numpy.zeros((centers.shape[0], 3))
    targets[:, :2] = centers

    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
    numBad = vi.findPoints(targets, tol2=1.e-10)
    assert numBad == 0, f'{numBad} cell-centre target points fell outside the grid'

    numeric = vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
    exact = exactGradient(centers[:, 0], centers[:, 1])
    return centers, numeric, exact


@pytest.mark.parametrize("M", [4, 8, 16, 32])
def test_vector_interp_reconstructs_gradient(M):
    """
    Part (3), per-resolution: the reconstructed vector should be close
    to the exact grad(chi) at every cell centre, with both the typical
    (RMS) and worst-case (max) error within a resolution-dependent
    tolerance -- see the module docstring for why these use different
    convergence orders (RMS ~ 1/M^2, max ~ 1/M, the latter set by the
    pole-adjacent cells).
    """
    centers, numeric, exact = _reconstructVectors(M)
    err = numpy.linalg.norm(numeric - exact, axis=-1)
    rms = numpy.sqrt((err ** 2).mean())

    tol_rms = 0.03 / M ** 2
    tol_max = 0.02 / M
    print(f'M={M}: max_err={err.max():.3e} (tol {tol_max:.3e}), '
          f'rms_err={rms:.3e} (tol {tol_rms:.3e})')
    assert rms < tol_rms, f'RMS error {rms:.3e} exceeds the O(1/M^2) tolerance {tol_rms:.3e}'
    assert err.max() < tol_max, f'max error {err.max():.3e} exceeds the O(1/M) tolerance {tol_max:.3e}'


def test_vector_interp_convergence():
    """
    Part (3), across resolutions: both the RMS and max error should
    shrink monotonically as the cubed-sphere is refined.
    """
    Ms = (4, 8, 16, 32)
    rms_errors, max_errors = [], []
    for M in Ms:
        centers, numeric, exact = _reconstructVectors(M)
        err = numpy.linalg.norm(numeric - exact, axis=-1)
        rms_errors.append(numpy.sqrt((err ** 2).mean()))
        max_errors.append(err.max())
        print(f'M={M}: max_err={max_errors[-1]:.3e}, rms_err={rms_errors[-1]:.3e}')

    for M_prev, M_next, e_prev, e_next in zip(Ms[:-1], Ms[1:], rms_errors[:-1], rms_errors[1:]):
        assert e_next < e_prev, (
            f'refining M={M_prev} -> M={M_next} should reduce the RMS error '
            f'({e_prev:.3e} -> {e_next:.3e}): the reconstruction is not converging')
    for M_prev, M_next, e_prev, e_next in zip(Ms[:-1], Ms[1:], max_errors[:-1], max_errors[1:]):
        assert e_next < e_prev, (
            f'refining M={M_prev} -> M={M_next} should reduce the max error '
            f'({e_prev:.3e} -> {e_next:.3e}): the reconstruction is not converging')


if __name__ == '__main__':
    test_exact_gradient_matches_finite_difference()
    for M in (4, 8, 16, 32):
        test_vector_interp_reconstructs_gradient(M)
    test_vector_interp_convergence()
