"""
Test mimetic 1-form interpolation for a velocity field that derives from a
vector potential (rigid-body rotation about an axis), on the sphere.

Field construction
-------------------
psi(r) = -Omega * (axis . r)     (a degree-1 solid spherical harmonic)
v      = rhat x grad(psi) = Omega * (axis x r)

psi and v are restrictions of *linear* functions of the ambient Cartesian
coordinates, hence C-infinity everywhere on the sphere, poles included --
that is what guarantees the regularity condition at the poles. AXIS is
deliberately tilted away from the grid's pole so that v is nonzero and
non-axisymmetric right where the (lon, lat) chart is singular, stress
testing the pole treatment (cf. the Williamson et al. 1992 solid-body
rotation advection test case).

Unlike the existing exact-potential / winding-form tests in this test
suite (test_polyline_integral.py, test_singularity2.py), this 1-form is
neither closed nor exact (d(v-flat) = 2*Omega*(axis . rhat) != 0), so it
exercises the genuinely rotational class of 1-forms rather than gradient
fields.

Exact line integral
--------------------
For r(t) tracing a great-circle arc from r0 to r1 at constant angular
speed, r x rdot = chat is conserved along the arc (angular momentum of
motion confined to a great circle), where

    chat = (r0 x r1) / |r0 x r1|

Hence v . rdot = Omega * axis . (r x rdot) = Omega * (axis . chat) is
itself constant along the arc, and integrates to the closed form

    integral_{r0}^{r1} v . dl = Omega * (axis . chat) * ds

where ds = angle between r0 and r1 (radians). This holds edge-by-edge
(used to build the exact discrete data below) and for whole great-circle
arcs used as test polylines, including ones that cross directly over a
pole.
"""
import numpy
from pathlib import Path
from mint import Grid, PolylineIntegral, NUM_EDGES_PER_QUAD, CELL_BY_CELL_DATA

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')
DEG2RAD = numpy.pi / 180.0

OMEGA = 1.0
# rotation axis tilted 45 degrees off the grid's pole
AXIS = numpy.array([numpy.sin(numpy.pi/4), 0.0, numpy.cos(numpy.pi/4)])

RESOLUTIONS = ('20x10', '40x20', '80x40')


def lonlat2xyz(lonlat):
    """
    Convert lon/lat (degrees, array with a trailing size-2 or size-3 axis)
    to unit Cartesian vectors.
    """
    lonlat = numpy.asarray(lonlat, dtype=numpy.float64)
    lon = lonlat[..., 0] * DEG2RAD
    lat = lonlat[..., 1] * DEG2RAD
    cosLat = numpy.cos(lat)
    xyz = numpy.empty(lon.shape + (3,), numpy.float64)
    xyz[..., 0] = cosLat * numpy.cos(lon)
    xyz[..., 1] = cosLat * numpy.sin(lon)
    xyz[..., 2] = numpy.sin(lat)
    return xyz


def arcIntegral(p0, p1, axis=AXIS, omega=OMEGA):
    """
    Exact integral of v . dl along the minor great-circle arc(s) from p0 to
    p1 (lon, lat in degrees; trailing axis of size 2, batchable over cells).
    """
    r0 = lonlat2xyz(p0)
    r1 = lonlat2xyz(p1)
    cross = numpy.cross(r0, r1)
    sinDs = numpy.linalg.norm(cross, axis=-1)
    cosDs = numpy.sum(r0 * r1, axis=-1)
    ds = numpy.arctan2(sinDs, cosDs)
    sinDsSafe = numpy.where(sinDs < 1.e-14, 1.0, sinDs)
    axisDotChat = numpy.dot(cross, axis) / sinDsSafe
    return omega * axisDotChat * ds


def slerp(p0, p1, n):
    """
    Sample n+1 points (lon, lat, 0) along the minor great-circle arc from
    p0 to p1 (lon, lat in degrees).
    """
    r0 = lonlat2xyz(p0)
    r1 = lonlat2xyz(p1)
    cosDs = numpy.dot(r0, r1)
    sinDs = numpy.linalg.norm(numpy.cross(r0, r1))
    ds = numpy.arctan2(sinDs, cosDs)
    t = numpy.linspace(0., 1., n + 1)
    if ds < 1.e-12:
        r = numpy.tile(r0, (n + 1, 1))
    else:
        a = numpy.sin((1. - t) * ds) / numpy.sin(ds)
        b = numpy.sin(t * ds) / numpy.sin(ds)
        r = numpy.outer(a, r0) + numpy.outer(b, r1)
        r /= numpy.linalg.norm(r, axis=1, keepdims=True)
    lon = numpy.arctan2(r[:, 1], r[:, 0]) / DEG2RAD
    lat = numpy.arctan2(r[:, 2], numpy.sqrt(r[:, 0]**2 + r[:, 1]**2)) / DEG2RAD
    xyz = numpy.zeros((n + 1, 3))
    xyz[:, 0] = lon
    xyz[:, 1] = lat
    return xyz


def geodesicPath(waypoints, nPerSeg):
    """
    Build a finely resolved polyline that follows the true great-circle
    geodesics through a sequence of (lon, lat) waypoints.
    """
    segs = [slerp(waypoints[i], waypoints[i + 1], nPerSeg)[:-1]
            for i in range(len(waypoints) - 1)]
    last = numpy.array([[waypoints[-1][0], waypoints[-1][1], 0.]])
    return numpy.vstack(segs + [last])


def buildEdgeData(grid):
    """
    Build the exact cell-by-cell edge-integrated data for v on the grid's
    own edges (used as the mimetic scheme's input field).
    """
    ncells = grid.getNumberOfCells()
    points = grid.getPoints()  # (ncells, 4, 3): lon, lat, elev in degrees
    data = numpy.zeros((ncells, NUM_EDGES_PER_QUAD))
    for i0 in range(NUM_EDGES_PER_QUAD):
        i1 = (i0 + 1) % NUM_EDGES_PER_QUAD
        # edges 2 and 3 run backwards relative to the (i0 -> i1) node order,
        # see mint.NUM_EDGES_PER_QUAD indexing convention
        sign = 1 - 2 * (i0 // 2)
        data[:, i0] = sign * arcIntegral(points[:, i0, :2], points[:, i1, :2])
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


def test_polar_crossing_arc():
    """
    A short great-circle arc that crosses directly over a pole of the
    (lon, lat) grid, checked against the exact closed-form arc integral.
    """
    p0 = (30.0, 85.0)
    p1 = (210.0, 85.0)  # antipodal in longitude -> minor arc runs through the near pole
    exact = arcIntegral(p0, p1)
    path = geodesicPath([p0, p1], nPerSeg=40)

    errors = []
    for res in RESOLUTIONS:
        flux = _computeFlux(res, path)
        error = abs(flux - exact)
        errors.append(error)
        print(f'polar_crossing_arc {res}: flux={flux:.6f} exact={exact:.6f} error={error:.2e}')

    # error should shrink as the mesh is refined
    assert errors[-1] < errors[0]


def test_pole_to_pole_loop():
    """
    Closed loop around the full meridian great circle through both poles
    (total angle swept = 2*pi), checked against the sum of the exact arc
    contributions. Routed through the equator at each meridian crossing to
    avoid ever connecting two antipodal points directly (undefined arc).
    """
    lon0 = 10.0
    waypoints = [(lon0, 90.0), (lon0, 0.0), (lon0, -90.0),
                 (lon0 + 180.0, 0.0), (lon0, 90.0)]
    exact = sum(arcIntegral(waypoints[i], waypoints[i + 1])
                for i in range(len(waypoints) - 1))
    path = geodesicPath(waypoints, nPerSeg=30)

    errors = []
    for res in RESOLUTIONS:
        flux = _computeFlux(res, path)
        error = abs(flux - exact)
        errors.append(error)
        print(f'pole_to_pole_loop {res}: flux={flux:.6f} exact={exact:.6f} error={error:.2e}')

    assert errors[-1] < errors[0]


if __name__ == '__main__':
    test_polar_crossing_arc()
    test_pole_to_pole_loop()
