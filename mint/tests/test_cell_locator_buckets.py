"""
Regression test for the lon-lat cell locator's bucket assignment. With a fixed
5 buckets in x, the number of buckets in y grows like
numCells/(5*numCellsPerBucket) and can exceed the number of cell rows. A face
used to be registered only in the buckets containing its nodes, so a point in
the interior of a face could land in a bucket listing none of the faces around
it: on a 360x180 lat-lon grid with numCellsPerBucket=10, no cell centroid was
found.
"""
import numpy
import pytest

from mint import Grid, VectorInterp


def create_latlon_points(nx, ny):
    lon = numpy.linspace(0., 360., nx + 1)
    lat = numpy.linspace(-90., 90., ny + 1)
    lon0, lat0 = numpy.meshgrid(lon[:-1], lat[:-1])
    lon1, lat1 = numpy.meshgrid(lon[1:], lat[1:])
    lon0, lat0, lon1, lat1 = (a.ravel() for a in (lon0, lat0, lon1, lat1))
    zero = numpy.zeros_like(lon0)
    points = numpy.empty((nx*ny, 4, 3), numpy.float64)
    points[:, 0, :] = numpy.stack([lon0, lat0, zero], axis=-1)
    points[:, 1, :] = numpy.stack([lon1, lat0, zero], axis=-1)
    points[:, 2, :] = numpy.stack([lon1, lat1, zero], axis=-1)
    points[:, 3, :] = numpy.stack([lon0, lat1, zero], axis=-1)
    return points


@pytest.mark.parametrize("nx, ny, numCellsPerBucket", [
    (360, 180, 10),   # far more buckets in y than cell rows
    (72, 36, 10),
    (100, 50, 10),
    (360, 180, 128),  # buckets larger than cells, worked before
])
def test_find_points_in_latlon_cells(nx, ny, numCellsPerBucket):

    points = create_latlon_points(nx, ny)

    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0)
    grid.setPoints(points)

    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=numCellsPerBucket, periodX=360.)

    # cell centroids
    centroids = points.mean(axis=1).copy()
    assert vi.findPoints(centroids, tol2=1.e-10) == 0

    # random points strictly inside the domain
    rng = numpy.random.default_rng(1234)
    npts = 2000
    targets = numpy.zeros((npts, 3), numpy.float64)
    targets[:, 0] = rng.uniform(0., 360., npts)
    targets[:, 1] = rng.uniform(-89.9, 89.9, npts)
    assert vi.findPoints(targets, tol2=1.e-10) == 0
