"""
Regression tests for numpy arrays handed over to the C library. C reads an
array as a dense block in row major (C) order, so arrays that are Fortran
ordered (e.g. a transpose) or not contiguous (e.g. a strided slice) used to be
read in the wrong order, silently giving wrong results. They must be rejected.
"""
import ctypes
from pathlib import Path

import numpy
import pytest

from mint import (Grid, PolylineIntegral, RegridEdges, VectorInterp, ExtensiveFieldAdaptor,
                  CELL_BY_CELL_DATA, UNIQUE_EDGE_DATA, FUNC_SPACE_W1, NUM_EDGES_PER_QUAD)

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def fortran_order(a):
    """Same values as a, Fortran ordered (not C contiguous unless 1D)"""
    return numpy.asfortranarray(a)


def strided(a):
    """Same values as a, as a non contiguous view into a larger array"""
    big = numpy.zeros(a.shape[:-1] + (2*a.shape[-1],), a.dtype)
    big[..., ::2] = a
    return big[..., ::2]


def every_other(a):
    """Same values as a 1D array, as a non contiguous view"""
    big = numpy.zeros(2*a.size, a.dtype)
    big[::2] = a.ravel()
    return big[::2]


def load_grid():
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/latlon4x2.nc$latlon')
    return grid


def test_helpers():
    a = numpy.arange(12.).reshape(4, 3)
    assert not fortran_order(a).flags.c_contiguous and (fortran_order(a) == a).all()
    assert not strided(a).flags.c_contiguous and (strided(a) == a).all()
    b = numpy.arange(5.)
    assert not every_other(b).flags.c_contiguous and (every_other(b) == b).all()


@pytest.mark.parametrize("make_bad", [fortran_order, strided])
def test_polyline_compute_weights(make_bad):
    grid = load_grid()
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=4, periodX=360.)
    xyz = numpy.array([(10., -20., 0.), (200., 30., 0.), (300., 60., 0.)])
    pli.computeWeights(xyz)
    with pytest.raises(ctypes.ArgumentError):
        pli.computeWeights(make_bad(xyz))


@pytest.mark.parametrize("make_bad", [fortran_order, strided])
def test_polyline_get_integral(make_bad):
    grid = load_grid()
    pli = PolylineIntegral()
    pli.setGrid(grid)
    pli.buildLocator(numCellsPerBucket=4, periodX=360.)
    pli.computeWeights(numpy.array([(10., -20., 0.), (200., 30., 0.)]))
    data = numpy.arange(grid.getNumberOfCells()*NUM_EDGES_PER_QUAD, dtype=numpy.float64)
    data = data.reshape(grid.getNumberOfCells(), NUM_EDGES_PER_QUAD)
    pli.getIntegral(data, CELL_BY_CELL_DATA)
    with pytest.raises(ctypes.ArgumentError):
        pli.getIntegral(make_bad(data), CELL_BY_CELL_DATA)


@pytest.mark.parametrize("which", ['src', 'dst'])
def test_regrid_apply(which):
    rg = RegridEdges()
    src, dst = load_grid(), load_grid()
    rg.setSrcGrid(src)
    rg.setDstGrid(dst)
    rg.buildLocator(numCellsPerBucket=4, periodX=360.)
    rg.computeWeights()
    src_data = numpy.arange(src.getNumberOfEdges(), dtype=numpy.float64)
    dst_data = numpy.zeros(dst.getNumberOfEdges())
    rg.apply(src_data, dst_data, placement=UNIQUE_EDGE_DATA)
    if which == 'src':
        src_data = every_other(src_data)
    else:
        dst_data = every_other(dst_data)
    with pytest.raises(ctypes.ArgumentError):
        rg.apply(src_data, dst_data, placement=UNIQUE_EDGE_DATA)


@pytest.mark.parametrize("make_bad", [fortran_order, strided])
def test_vectorinterp_find_points(make_bad):
    grid = load_grid()
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=4, periodX=360.)
    targets = numpy.array([(10., -20., 0.), (200., 30., 0.), (300., 60., 0.)])
    assert vi.findPoints(targets, tol2=1.e-10) == 0
    with pytest.raises(ctypes.ArgumentError):
        vi.findPoints(make_bad(targets), tol2=1.e-10)


@pytest.mark.parametrize("make_bad", [fortran_order, strided])
def test_vectorinterp_get_edge_vectors(make_bad):
    grid = load_grid()
    vi = VectorInterp()
    vi.setGrid(grid)
    vi.buildLocator(numCellsPerBucket=4, periodX=360.)
    vi.findPoints(numpy.array([(10., -20., 0.), (200., 30., 0.)]), tol2=1.e-10)
    data = numpy.ones((grid.getNumberOfCells(), NUM_EDGES_PER_QUAD))
    vi.getEdgeVectors(data, placement=CELL_BY_CELL_DATA)
    with pytest.raises(ctypes.ArgumentError):
        vi.getEdgeVectors(make_bad(data), placement=CELL_BY_CELL_DATA)


def test_extensive_field_adaptor():
    grid = load_grid()
    efa = ExtensiveFieldAdaptor()
    efa.setGrid(grid)
    n = grid.getNumberOfEdges()
    u = numpy.ones(n)
    v = numpy.ones(n)
    data = numpy.zeros(n)
    efa.fromVectorField(u, v, data, UNIQUE_EDGE_DATA, FUNC_SPACE_W1)
    with pytest.raises(ctypes.ArgumentError):
        efa.fromVectorField(every_other(u), v, data, UNIQUE_EDGE_DATA, FUNC_SPACE_W1)


@pytest.mark.parametrize("make_bad", [fortran_order, strided])
def test_grid_set_points(make_bad):
    points = load_grid().getPoints().copy()
    grid = Grid()
    grid.setPoints(points)
    with pytest.raises((ctypes.ArgumentError, RuntimeError)):
        Grid().setPoints(make_bad(points))
