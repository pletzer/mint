"""
Line integrals and edge regridding with useXYZLocator=True on a whole sphere:
the C24 FV3 cubed-sphere grid given in unit-sphere (x, y, z) coordinates.

For edge data that are exact differences of a potential, MINT's results must
telescope across cell boundaries: closed paths and regridded target cells have
zero circulation to round-off, whatever the path and including at the poles,
across panel edges and along source edges (lon -10 and 170 are panel
boundaries of the FV3 grid). Every target cell must also receive data.
"""
from pathlib import Path

import numpy
import pytest

from mint import Grid, PolylineIntegral, RegridEdges, CELL_BY_CELL_DATA

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

# (start, end) node of each cell-by-cell edge
EDGE_NODES = ((0, 1), (1, 2), (3, 2), (0, 3))


def lonlat_to_xyz(lonlat):
    lon, lat = numpy.deg2rad(lonlat[..., 0]), numpy.deg2rad(lonlat[..., 1])
    return numpy.stack([numpy.cos(lat)*numpy.cos(lon),
                        numpy.cos(lat)*numpy.sin(lon),
                        numpy.sin(lat)], axis=-1)


def potential(xyz):
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    return x + x*z + 0.5*y


def edge_data(points):
    phi = potential(points)
    return numpy.ascontiguousarray(
        numpy.stack([phi[:, b] - phi[:, a] for a, b in EDGE_NODES], axis=1).ravel())


def xyz_grid(points):
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=False)
    grid.setPoints(points)
    return grid


@pytest.fixture(scope='module')
def c24():
    g = Grid()
    g.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1, degrees=True)
    g.loadFromFV32DFile(str(DATA_DIR / 'C24_SCRIP_desc.181018.nc'))
    points = numpy.ascontiguousarray(lonlat_to_xyz(g.getPoints()))
    return points, xyz_grid(points), edge_data(points)


@pytest.mark.parametrize("numCellsPerBucket", [64, 129])
@pytest.mark.parametrize("centre,radius,nseg", [
    ((-10., 5.), 1., 8),        # inside a panel
    ((-10., 5.), 15., 32),      # across several cells
    ((-55., 0.), 5., 32),       # across a panel edge
    ((35., 89.), 3., 16),       # around the north pole
    ((-10., -60.), 3., 16),
])
def test_closed_path(c24, numCellsPerBucket, centre, radius, nseg):
    _, grid, data = c24
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=numCellsPerBucket, useXYZLocator=True)
    g = numpy.linspace(0., 2*numpy.pi, nseg + 1)
    lonlat = numpy.stack([centre[0] + radius*numpy.sin(g),
                          centre[1] + radius*numpy.cos(g)], axis=-1)
    lonlat[-1] = lonlat[0]
    pli.computeWeights(numpy.ascontiguousarray(lonlat_to_xyz(lonlat)), counterclock=False)
    circulation = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)
    print(f'circulation = {circulation}')
    assert abs(circulation) < 1.e-12


def test_open_path(c24):
    _, grid, data = c24
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=128, useXYZLocator=True)
    # panel centres, which are grid nodes: the result is exact there
    lonlat = numpy.array([[-100., 0.], [-10., 0.]])
    xyz = numpy.ascontiguousarray(lonlat_to_xyz(lonlat))
    pli.computeWeights(xyz, counterclock=False)
    integral = pli.getIntegral(data=data, placement=CELL_BY_CELL_DATA)
    exact = potential(xyz[1]) - potential(xyz[0])
    print(f'integral = {integral} exact = {exact}')
    assert abs(integral - exact) < 1.e-12


def test_regrid_to_lonlat(c24):
    _, srcGrid, srcData = c24
    # 10 x 10 deg lon-lat target, including the degenerate cells at the poles
    lon = numpy.linspace(0., 360., 37)
    lat = numpy.linspace(-90., 90., 19)
    LON0, LAT0 = numpy.meshgrid(lon[:-1], lat[:-1], indexing='ij')
    LON1, LAT1 = numpy.meshgrid(lon[1:], lat[1:], indexing='ij')
    lonlat = numpy.stack([numpy.stack([LON0, LAT0], axis=-1),
                          numpy.stack([LON1, LAT0], axis=-1),
                          numpy.stack([LON1, LAT1], axis=-1),
                          numpy.stack([LON0, LAT1], axis=-1)], axis=2).reshape(-1, 4, 2)
    dstPoints = numpy.ascontiguousarray(lonlat_to_xyz(lonlat))
    dstGrid = xyz_grid(dstPoints)

    rg = RegridEdges()
    rg.setSrcGrid(srcGrid)
    rg.setDstGrid(dstGrid)
    rg.buildLocator(numCellsPerBucket=128, useXYZLocator=True)
    rg.computeWeights()
    dstData = numpy.zeros(4*dstPoints.shape[0])
    rg.apply(srcData, dstData, placement=CELL_BY_CELL_DATA)

    e = dstData.reshape(-1, 4)
    circulation = e[:, 0] + e[:, 1] - e[:, 2] - e[:, 3]
    exact = edge_data(dstPoints).reshape(-1, 4)
    missing = numpy.all(e == 0, axis=1) & numpy.any(exact != 0, axis=1)
    err = numpy.abs(e - exact).max()
    print(f'max circulation = {numpy.abs(circulation).max()}, max edge error = {err}, '
          f'missing cells = {missing.sum()}')
    assert numpy.abs(circulation).max() < 1.e-12
    assert missing.sum() == 0
    # end points of each target edge are grid nodes or inside cells: O(h^2)
    assert err < 1.e-2
