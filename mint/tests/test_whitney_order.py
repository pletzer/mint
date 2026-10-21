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


def structuredGrid(nx, ny):
    """flat nx x ny UGRID grid with nodes at the integer points (i, j)"""
    node = lambda i, j: j*(nx + 1) + i
    xyz = numpy.array([(i, j, 0.) for j in range(ny + 1) for i in range(nx + 1)], numpy.float64)
    f2n = numpy.array([(node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1))
                       for j in range(ny) for i in range(nx)], numpy.uint64)
    e2n = numpy.array([(node(i, j), node(i + 1, j)) for j in range(ny + 1) for i in range(nx)] +
                      [(node(i, j), node(i, j + 1)) for j in range(ny) for i in range(nx + 1)],
                      numpy.uint64)
    grid = Grid()
    grid.setFlags(0, 0, degrees=False)
    grid.loadFromUgrid2DData(xyz, f2n, e2n)
    return grid, xyz, f2n


def test_whitney_element():
    """
    Covariant components from the element: the gradient of a quadratic potential
    is reproduced exactly by order 2, not by order 1
    """
    nx, ny = 6, 4
    grid, xyz, f2n = structuredGrid(nx, ny)

    def phi(s, t):
        return s*s + 0.5*s*t - t*t

    # cell by cell edge integrals, edges in the positive xi/eta direction
    pts = xyz[f2n.astype(int)]                      # (ncells, 4, 3)
    data = numpy.zeros((nx*ny, NUM_EDGES_PER_QUAD))
    for ie, (i0, i1) in enumerate(EDGE_VERTS):
        data[:, ie] = phi(pts[:, i1, 0], pts[:, i1, 1]) - phi(pts[:, i0, 0], pts[:, i0, 1])
    data = data.ravel()

    rng = numpy.random.default_rng(7)
    cellIds = rng.integers(0, nx*ny, 50)
    xi, eta = rng.random(50), rng.random(50)
    s = pts[cellIds, 0, 0] + xi
    t = pts[cellIds, 0, 1] + eta
    dphids, dphidt = 2*s + 0.5*t, 0.5*s - 2*t

    e1 = mint.WhitneyElement(1, grid)
    assert e1.getOrder() == 1
    a, b = e1.getCovariantComponents(cellIds, xi, eta, data)
    # lowest order: a = (1 - eta) d0 + eta d2, b = (1 - xi) d3 + xi d1
    d = data.reshape(-1, 4)[cellIds]
    assert numpy.allclose(a, (1 - eta)*d[:, 0] + eta*d[:, 2], atol=1e-14)
    assert numpy.allclose(b, (1 - xi)*d[:, 3] + xi*d[:, 1], atol=1e-14)
    assert numpy.abs(a - dphids).max() > 0.1

    for biquadratic in (True, False):
        e2 = mint.WhitneyElement(2, grid, biquadratic=biquadratic)
        assert e2.getOrder() == 2
        a, b = e2.getCovariantComponents(cellIds, xi, eta, data)
        assert numpy.abs(a - dphids).max() < 1e-12
        assert numpy.abs(b - dphidt).max() < 1e-12

    # invalid cells are reported
    with pytest.raises(RuntimeError):
        e2.getCovariantComponents([nx*ny], [0.5], [0.5], data)
    # order 2 needs a grid that can be tiled into 2x2 blocks
    oddGrid = structuredGrid(5, 4)[0]
    with pytest.raises(RuntimeError):
        mint.WhitneyElement(2, oddGrid)
