"""
Line integrals with PolylineIntegral.buildLocator(useXYZLocator=True), on a flat
grid lying in a plane of 3D space: a tilted plane and a vertical (x = const) one,
whose projection onto (x, y) is a line. The integral of grad(phi) along a polyline
of a linear potential phi is phi(end) - phi(start), provided the segments traced
through the cells cover the whole polyline.
"""
from mint import Grid, PolylineIntegral, CELL_BY_CELL_DATA
import numpy
import pytest


ORIGIN = numpy.array([1., 1., 1.])

# orthonormal bases (e1, e2) of the planes
PLANES = {
    # normal (1, 1, 1)/sqrt(3)
    'tilted': (numpy.array([1., -1., 0.]) / numpy.sqrt(2.),
               numpy.cross(numpy.array([1., 1., 1.]) / numpy.sqrt(3.),
                           numpy.array([1., -1., 0.]) / numpy.sqrt(2.))),
    # normal (1, 0, 0)
    'x_const': (numpy.array([0., 1., 0.]), numpy.array([0., 0., 1.])),
}


def to_xyz(u, v, plane='tilted'):
    e1, e2 = PLANES[plane]
    return ORIGIN + u*e1 + v*e2


def potential(p):
    x, y, z = p
    return x + 2.*y - 3.*z


def build_grid_and_data(nx, ny, plane='tilted'):
    points = numpy.zeros((nx*ny, 4, 3), numpy.float64)
    data = numpy.zeros((nx*ny, 4), numpy.float64)
    k = 0
    for j in range(ny):
        for i in range(nx):
            u0, u1 = i/nx, (i + 1)/nx
            v0, v1 = j/ny, (j + 1)/ny
            # 3-->--2
            # |     |
            # ^     ^
            # |     |
            # 0-->--1
            points[k, 0, :] = to_xyz(u0, v0, plane)
            points[k, 1, :] = to_xyz(u1, v0, plane)
            points[k, 2, :] = to_xyz(u1, v1, plane)
            points[k, 3, :] = to_xyz(u0, v1, plane)
            phi = [potential(points[k, i, :]) for i in range(4)]
            data[k, :] = phi[1] - phi[0], phi[2] - phi[1], phi[2] - phi[3], phi[3] - phi[0]
            k += 1
    return points, data


@pytest.mark.parametrize("plane", list(PLANES))
@pytest.mark.parametrize("uv", [
    [(0.05, 0.03), (0.95, 0.87)],
    [(0.05, 0.55), (0.90005, 0.55)],               # ends just past a cell boundary
    [(0.1, 0.1), (0.9, 0.2), (0.5, 0.9), (0.12, 0.6)],
])
def test_polyline_integral_xyz_plane(plane, uv):

    points, data = build_grid_and_data(10, 10, plane)

    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=False)
    grid.setPoints(points)

    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=4, useXYZLocator=True)

    xyz = numpy.array([to_xyz(u, v, plane) for u, v in uv])
    pli.computeWeights(xyz, counterclock=False)
    flux = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)

    exact = potential(xyz[-1]) - potential(xyz[0])
    print(f'flux = {flux} exact = {exact} error = {flux - exact}')
    assert abs(flux - exact) < 1.e-10
