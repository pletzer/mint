"""
Test mimetic 1-form interpolation for a velocity field that derives from a
vector potential engineered to be finite and nonzero at one pole and exactly
zero at the other.

Field construction
-------------------
psi(lambda, theta) = cos(theta) * (1 + sin(theta)) * cos(lambda)
                    = x * (1 + z)     (Cartesian x = cos(theta)cos(lambda),
                                       z = sin(theta), on the unit sphere)

u = curl(psi * rhat) has physical components (unit sphere; this line
integral is independent of the sphere radius, see below)

    u_lambda = (cos(2*theta) - sin(theta)) * cos(lambda)
    u_theta  = (1 + sin(theta)) * sin(lambda)

psi is an azimuthal wavenumber-1 (m = 1) mode: only |m| = 1 modes can have a
finite, nonzero velocity right at a pole (m = 0 always vanishes there, |m|
>= 2 always vanishes there too). Its envelope g(sin(theta)) = 1 + sin(theta)
is chosen so g(1) = 2 != 0 = g(-1): the north pole has finite nonzero
velocity (0, -2, 0), the south pole has zero velocity -- something a zonal
or single tilted-axis (rigid rotation) field can never do, since a rotation
axis's two endpoints are always antipodal, forcing both geographic poles to
be zero together or nonzero together.

Unlike test_pole_rotation.py's rigid-body rotation field, whose line
integral along a great-circle arc has a clean closed form via conservation
of r x rdot, this field has no comparably simple invariant for an arbitrary
path. Its line integral IS however exactly computable along a path that is
straight in (lambda, theta) -- see the "Regularity of the velocity field at
the north pole" write-up accompanying work/vector_potential -- and that
exact value is reproduced here by a fine quadrature (straightLineIntegral),
used both as the ground truth for the test polylines and, edge by edge, to
build the mimetic scheme's exact input data (buildEdgeData), exactly as
test_pole_rotation.py does with its own (approximate, great-circle-based)
arcIntegral.
"""
import numpy
from pathlib import Path
from mint import Grid, PolylineIntegral, NUM_EDGES_PER_QUAD, CELL_BY_CELL_DATA

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
DEG2RAD = numpy.pi / 180.0

RESOLUTIONS = ('20x10', '40x20', '80x40')


def uDotDl(lam, theta, dlam, dtheta):
    """
    u . dl for the pole-asymmetric field, at points (lam, theta) along a
    path with local displacement (dlam, dtheta) -- i.e. the integrand of
    the line integral before integrating over the path parameter.
    """
    return (numpy.cos(theta) * (numpy.cos(2. * theta) - numpy.sin(theta)) * numpy.cos(lam) * dlam
            + (1. + numpy.sin(theta)) * numpy.sin(lam) * dtheta)


def straightLineIntegral(p0, p1, n=50):
    """
    Line integral of u . dl along the path STRAIGHT IN (lon, lat) from p0
    to p1 (lon, lat in degrees; trailing axis of size 2, batchable over
    cells), evaluated by an n-point quadrature. Exact up to quadrature
    error, which is negligible for the short segments (grid edges) this is
    used on; a much finer quadrature is used separately below for the
    ground-truth value of a whole test polyline.
    """
    p0 = numpy.asarray(p0, dtype=numpy.float64)
    p1 = numpy.asarray(p1, dtype=numpy.float64)
    lam0, th0 = p0[..., 0] * DEG2RAD, p0[..., 1] * DEG2RAD
    lam1, th1 = p1[..., 0] * DEG2RAD, p1[..., 1] * DEG2RAD
    dlam, dth = lam1 - lam0, th1 - th0

    tt = numpy.linspace(0., 1., n + 1)
    lam = lam0[..., None] + tt * dlam[..., None]
    theta = th0[..., None] + tt * dth[..., None]
    integrand = uDotDl(lam, theta, dlam[..., None], dth[..., None])

    dt = tt[1] - tt[0]
    return dt * (integrand[..., 0] / 2. + integrand[..., 1:-1].sum(axis=-1) + integrand[..., -1] / 2.)


def buildEdgeData(grid):
    """
    Build the (quadrature-)exact cell-by-cell edge-integrated data for u on
    the grid's own edges (used as the mimetic scheme's input field). Grid
    edges are straight in (lon, lat), so straightLineIntegral is exact
    (to quadrature accuracy) for them, unlike test_pole_rotation.py's
    great-circle-based arcIntegral applied to the same kind of edges.
    """
    ncells = grid.getNumberOfCells()
    points = grid.getPoints()  # (ncells, 4, 3): lon, lat, elev in degrees
    data = numpy.zeros((ncells, NUM_EDGES_PER_QUAD))
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        # edges 2 and 3 run backwards relative to the (i0 -> i1) node order,
        # see mint.NUM_EDGES_PER_QUAD indexing convention
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * straightLineIntegral(points[:, i0, :2], points[:, i1, :2])
    return data


def _computeFlux(res, xyzPath):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/latlon{res}Shifted.nc$mesh')
    data = buildEdgeData(grid)

    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
    pli.computeWeights(xyzPath, counterclock=False)
    return pli.getIntegral(data, placement=CELL_BY_CELL_DATA)


def test_line_integral_to_north_pole():
    """
    A path straight in (lon, lat) ending exactly at the north pole, where
    the field is finite and nonzero. Checked against the quadrature-exact
    straight-line integral, at increasing grid resolution.
    """
    p0 = (-30.0, -20.0)
    p1 = (45.0, 90.0)
    exact = straightLineIntegral(p0, p1, n=200000)
    path = numpy.array([(p0[0], p0[1], 0.), (p1[0], p1[1], 0.)])

    errors = []
    for res in RESOLUTIONS:
        flux = _computeFlux(res, path)
        error = abs(flux - exact)
        errors.append(error)
        print(f'to_north_pole {res}: flux={flux:.6f} exact={exact:.6f} error={error:.2e}')

    assert abs(exact) > 0.1, "sanity check: the field should be clearly nonzero at the north pole"
    # error should shrink as the mesh is refined
    assert errors[-1] < errors[0]


def test_line_integral_to_south_pole():
    """
    A path straight in (lon, lat) ending exactly at the south pole. Note
    this line integral itself is generally nonzero even though the FIELD
    is exactly zero right at the south pole endpoint (see
    test_near_pole_flux_asymmetry below for that): the integral
    accumulates u . dl along the whole path, most of which is away from
    the pole. Checked against the quadrature-exact straight-line integral,
    at increasing grid resolution.
    """
    p0 = (-30.0, -20.0)
    p1 = (45.0, -90.0)
    exact = straightLineIntegral(p0, p1, n=200000)
    path = numpy.array([(p0[0], p0[1], 0.), (p1[0], p1[1], 0.)])

    errors = []
    for res in RESOLUTIONS:
        flux = _computeFlux(res, path)
        error = abs(flux - exact)
        errors.append(error)
        print(f'to_south_pole {res}: flux={flux:.6f} exact={exact:.6f} error={error:.2e}')

    # error should shrink as the mesh is refined
    assert errors[-1] < errors[0]


def test_pole_values_are_asymmetric():
    """
    Direct (analytic, not MINT) check of the property motivating this test
    module: sampling u (via its line-integral effect over a tiny
    displacement) at the two poles shows the north pole is finite and
    nonzero while the south pole is exactly zero -- confirmed analytically
    via g(1) = 2 != 0 = g(-1) for g(s) = 1 + s.
    """
    north = uDotDl(lam=numpy.array([0., numpy.pi / 3, numpy.pi]),
                    theta=numpy.pi / 2, dlam=0., dtheta=1.)
    south = uDotDl(lam=numpy.array([0., numpy.pi / 3, numpy.pi]),
                    theta=-numpy.pi / 2, dlam=0., dtheta=1.)
    # u_theta(north pole) = (1 + sin(pi/2)) * sin(lambda) = 2*sin(lambda):
    # nonzero for generic lambda, exactly matching g(1) = 2
    assert numpy.any(numpy.abs(north) > 1.0)
    # u_theta(south pole) = (1 + sin(-pi/2)) * sin(lambda) = 0 identically
    assert numpy.allclose(south, 0.0)


def test_near_pole_flux_asymmetry():
    """
    The MINT-computed (not just analytic) counterpart of
    test_pole_values_are_asymmetric: a short meridional segment ending at
    each pole, at the longitude (90 deg) where the m=1 envelope is largest.
    Because u_theta ~ g(1) = 2 (finite) near the north pole but u_theta ~
    (eps_rad)^2/2 -> 0 near the south pole, a short segment's flux should
    be many orders of magnitude larger near the north pole than near the
    south pole, for the same angular length -- checked here through an
    actual PolylineIntegral computation, not just the formula.
    """
    eps = 5.0  # degrees
    pathNorth = numpy.array([(90.0, 90.0 - eps, 0.), (90.0, 90.0, 0.)])
    pathSouth = numpy.array([(90.0, -90.0 + eps, 0.), (90.0, -90.0, 0.)])

    exactNorth = straightLineIntegral((90.0, 90.0 - eps), (90.0, 90.0), n=200000)
    exactSouth = straightLineIntegral((90.0, -90.0 + eps), (90.0, -90.0), n=200000)
    assert abs(exactNorth) > 0.1
    assert abs(exactSouth) < 1.e-3

    for res in RESOLUTIONS:
        fluxNorth = _computeFlux(res, pathNorth)
        fluxSouth = _computeFlux(res, pathSouth)
        ratio = abs(fluxNorth / fluxSouth)
        print(f'near_pole_flux_asymmetry {res}: north={fluxNorth:.6f} (exact {exactNorth:.6f}), '
              f'south={fluxSouth:.3e} (exact {exactSouth:.3e}), |north/south|={ratio:.1f}')
        # the near-north flux should always dwarf the near-south flux for
        # the same short angular segment
        assert ratio > 100.


if __name__ == '__main__':
    test_line_integral_to_north_pole()
    test_line_integral_to_south_pole()
    test_pole_values_are_asymmetric()
    test_near_pole_flux_asymmetry()
