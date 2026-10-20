"""
Selecting the order of the edge element: 1 (lowest order Whitney element, the
default) or 2 (second order element on 2x2 blocks of grid cells, which uses the
same edge data).
"""
from mint import (Grid, VectorInterp, PolylineIntegral, RegridEdges,
                  UNIQUE_EDGE_DATA, NUM_EDGES_PER_QUAD)
import mint
import numpy
import pytest
from pathlib import Path
from iris_utils import _u_v_cubes_from_ugrid_file


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

# local vertex pairs of the cell edges, pointing in the positive xi/eta direction
EDGE_VERTS = [(0, 1), (1, 2), (3, 2), (0, 3)]


def loadGrid(fileAndMesh, flags):
    grid = Grid()
    grid.setFlags(*flags)
    grid.loadFromUgrid2DFile(str(DATA_DIR / fileAndMesh))
    return grid


def potential(lon, lat):
    """smooth potential on the sphere, lon and lat in degrees"""
    lam, the = numpy.radians(lon), numpy.radians(lat)
    x, y, z = numpy.cos(the)*numpy.cos(lam), numpy.cos(the)*numpy.sin(lam), numpy.sin(the)
    return x*y + 0.5*z*z - 0.3*x*z + 0.7*y + 0.2*numpy.sin(2*x + z)


def gradientEdgeData(grid):
    """exact unique edge integrals of the gradient of the potential"""
    pts = grid.getPoints()
    data = numpy.zeros(grid.getNumberOfEdges(), numpy.float64)
    for icell in range(grid.getNumberOfCells()):
        for ie, (i0, i1) in enumerate(EDGE_VERTS):
            edgeId, edgeSign = grid.getEdgeId(icell, ie)
            p0, p1 = pts[icell, i0, :], pts[icell, i1, :]
            data[edgeId] = edgeSign*(potential(p1[0], p1[1]) - potential(p0[0], p0[1]))
    return data


def test_set_order():

    grid = loadGrid('cs_16.nc$physics', (1, 1))

    vi = VectorInterp()
    # order 2 needs a grid
    with pytest.raises(RuntimeError):
        vi.setOrder(2)
    vi.setGrid(grid)
    assert vi.getOrder() == 1
    vi.setOrder(2)
    assert vi.getOrder() == 2
    # invalid orders are rejected and leave the order unchanged
    with pytest.raises(RuntimeError):
        vi.setOrder(3)
    assert vi.getOrder() == 2
    vi.setOrder(1)
    assert vi.getOrder() == 1

    pli = PolylineIntegral()
    pli.setGrid(grid)
    assert pli.getOrder() == 1
    pli.setOrder(2)
    assert pli.getOrder() == 2

    rg = RegridEdges()
    rg.setSrcGrid(grid)
    assert rg.getOrder() == 1
    rg.setOrder(2)
    assert rg.getOrder() == 2


def test_polyline_integral_order2():
    """
    The line integral of a gradient field is the difference of the discrete
    potential at the end points: O(h^2) with order 1, O(h^3) with order 2
    """
    grid = loadGrid('cs_64.nc$physics', (1, 1))
    data = gradientEdgeData(grid)

    rng = numpy.random.default_rng(3)
    paths = [numpy.stack([rng.uniform(-150., 150., 4), rng.uniform(-70., 70., 4),
                          numpy.zeros(4)], -1) for _ in range(20)]

    errors = {}
    for order in (1, 2):
        pli = PolylineIntegral()
        pli.setGrid(grid)
        pli.setOrder(order)
        pli.buildLocator(numCellsPerBucket=100, periodX=360., enableFolding=True)
        err = []
        for xyz in paths:
            pli.computeWeights(xyz)
            result = pli.getIntegral(data, placement=UNIQUE_EDGE_DATA)
            exact = potential(xyz[-1, 0], xyz[-1, 1]) - potential(xyz[0, 0], xyz[0, 1])
            err.append(abs(result - exact))
        errors[order] = numpy.sqrt(numpy.mean(numpy.square(err)))
    print(f'RMS line integral errors on C64: order 1 {errors[1]:.3e} order 2 {errors[2]:.3e}')
    assert errors[2] < 0.1*errors[1]


def test_regrid_order2():
    """
    Regridding the edge integrals of a gradient field from a lon-lat grid onto a
    cubed-sphere grid
    """
    src = loadGrid('latlon100x50.nc$latlon', (0, 0))
    dst = loadGrid('cs_16.nc$physics', (1, 1))
    srcData = gradientEdgeData(src)
    exact = gradientEdgeData(dst)

    errors = {}
    for order in (1, 2):
        rg = RegridEdges()
        rg.setSrcGrid(src)
        rg.setDstGrid(dst)
        rg.setOrder(order)
        rg.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
        rg.computeWeights()
        dstData = numpy.zeros(rg.getNumDstEdges(), numpy.float64)
        rg.apply(srcData, dstData, placement=UNIQUE_EDGE_DATA)
        errors[order] = numpy.sqrt(numpy.mean(numpy.square(dstData - exact)))
    print(f'RMS regridded edge integral errors: order 1 {errors[1]:.3e} order 2 {errors[2]:.3e}')
    assert errors[2] < 0.25*errors[1]


def test_iris_order():
    u, v = _u_v_cubes_from_ugrid_file(DATA_DIR / 'cs128_wind_zt.nc')
    flx = mint.IrisMintFlux(u.mesh, src_flags=(1, 1, 1), tgt_line=[(0., 0.), (90., 30.)], order=2)
    assert flx.flux_calc.getOrder() == 2
    flx = mint.IrisMintFlux(u.mesh, src_flags=(1, 1, 1), tgt_line=[(0., 0.), (90., 30.)])
    assert flx.flux_calc.getOrder() == 1
