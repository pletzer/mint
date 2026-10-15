"""
Regression tests for Python objects whose memory is borrowed by the C library:
grids set on RegridEdges/VectorInterp/PolylineIntegral/ExtensiveFieldAdaptor/
cell locators, locators set on VectorInterp, arrays attached to a Grid without
copying and the point array view returned by Grid.getPoints. The wrapper that
borrows the object must keep it alive -- without copying it -- for as long as
the C side may use it, and release it afterwards.
"""
import gc
import weakref
from pathlib import Path

import numpy

from mint import (Grid, RegridEdges, VectorInterp, PolylineIntegral, ExtensiveFieldAdaptor,
                  LonLatCellLocator, XYZCellLocator)

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def load_grid():
    grid = Grid()
    grid.setFlags(fixLonAcrossDateline=0, averageLonAtPole=0, degrees=True)
    grid.loadFromUgrid2DFile(f'{DATA_DIR}/latlon4x2.nc$latlon')
    return grid


def is_alive(ref):
    gc.collect()
    return ref() is not None


def test_regrid_edges_grids():
    rg = RegridEdges()
    src, dst = load_grid(), load_grid()
    wsrc, wdst = weakref.ref(src), weakref.ref(dst)
    rg.setSrcGrid(src)
    rg.setDstGrid(dst)
    del src, dst
    assert is_alive(wsrc) and is_alive(wdst)
    # the grids are still usable through the regridder
    rg.buildLocator(numCellsPerBucket=4, periodX=360.)
    rg.computeWeights()
    # replacing a grid releases the previous one
    rg.setSrcGrid(load_grid())
    assert not is_alive(wsrc) and is_alive(wdst)
    del rg
    assert not is_alive(wdst)


def test_vector_interp_grid_and_locator():
    vi = VectorInterp()
    grid = load_grid()
    wgrid = weakref.ref(grid)
    vi.setGrid(grid)
    loc = LonLatCellLocator()
    wloc = weakref.ref(loc)
    loc.setDataSet(grid)
    loc.setNumberOfCellsPerBucket(4)
    loc.buildLocator()
    vi.setLocator(loc)
    del grid, loc
    assert is_alive(wgrid) and is_alive(wloc)
    assert vi.findPoints(numpy.array([[45., 10., 0.]]), tol2=1.e-10) == 0
    # buildLocator replaces the borrowed locator with one owned by the C side
    vi.buildLocator(numCellsPerBucket=4, periodX=360.)
    assert not is_alive(wloc) and is_alive(wgrid)
    del vi
    assert not is_alive(wgrid)


def test_cell_locator_data_set():
    for cls in (LonLatCellLocator, XYZCellLocator):
        loc = cls()
        grid = load_grid()
        wgrid = weakref.ref(grid)
        loc.setDataSet(grid)
        del grid
        assert is_alive(wgrid)
        loc.setNumberOfCellsPerBucket(4)
        loc.buildLocator()
        del loc
        assert not is_alive(wgrid)


def test_polyline_integral_grid():
    pli = PolylineIntegral()
    grid = load_grid()
    wgrid = weakref.ref(grid)
    pli.setGrid(grid)
    del grid
    assert is_alive(wgrid)
    pli.buildLocator(numCellsPerBucket=4, periodX=360.)
    pli.computeWeights(numpy.array([(10., -20., 0.), (200., 30., 0.)]))
    del pli
    assert not is_alive(wgrid)


def test_extensive_field_adaptor_grid():
    efa = ExtensiveFieldAdaptor()
    grid = load_grid()
    wgrid = weakref.ref(grid)
    efa.setGrid(grid)
    del grid
    assert is_alive(wgrid)
    del efa
    assert not is_alive(wgrid)


def test_attach_without_copy():
    grid = load_grid()
    data = numpy.arange(grid.getNumberOfCells(), dtype=numpy.float64)
    wdata = weakref.ref(data)
    grid.attach('cellIds', data, copy=False)
    del data
    assert is_alive(wdata)
    del grid
    assert not is_alive(wdata)


def test_get_points_view():
    grid = load_grid()
    wgrid = weakref.ref(grid)
    points = grid.getPoints()
    expected = points.copy()
    del grid
    assert is_alive(wgrid)
    # still a view (no copy), still reading valid memory
    assert not points.flags.owndata
    assert (points == expected).all()
    del points
    assert not is_alive(wgrid)
