from mint import RegridEdges, Grid, UNIQUE_EDGE_DATA
import netCDF4
import numpy
import pytest
from pathlib import Path
from tempfile import TemporaryDirectory


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def createRegridder():
    src_grid = Grid()
    src_grid.setFlags(0, 0) # lon-lat
    src_grid.loadFromUgrid2DFile(f'{DATA_DIR / "latlon100x50.nc"}$latlon')
    dst_grid = Grid()
    dst_grid.setFlags(1, 1) # cubed-sphere
    dst_grid.loadFromUgrid2DFile(f'{DATA_DIR / "cs_16.nc"}$physics')
    rg = RegridEdges()
    rg.setSrcGrid(src_grid)
    rg.setDstGrid(dst_grid)
    # keep the grids alive as long as the regridder
    return rg, src_grid, dst_grid


def regrid(rg, src_data):
    dst_data = numpy.zeros(rg.getNumDstEdges(), numpy.float64)
    rg.apply(src_data, dst_data, placement=UNIQUE_EDGE_DATA)
    return dst_data


def test_weights_order_attributes():

    rg, sg, dg = createRegridder()
    rg.buildLocator(numCellsPerBucket=128, periodX=360., enableFolding=False)
    rg.computeWeights()

    src_data = numpy.random.default_rng(1).standard_normal(rg.getNumSrcEdges())
    expected = regrid(rg, src_data)

    with TemporaryDirectory() as d:

        filename = str(Path(d) / 'weights.nc')
        rg.dumpWeights(filename)
        with netCDF4.Dataset(filename) as nc:
            assert nc.src_element_order == 1
            assert nc.src_num_dofs_per_cell == 4
            # the edge description is unchanged
            assert nc.dimensions['num_edges_per_cell'].size == 4

        # reading back the weights gives the same result
        rg2, sg2, dg2 = createRegridder()
        rg2.loadWeights(filename)
        assert numpy.array_equal(regrid(rg2, src_data), expected)

        # files written before the attributes existed are lowest order
        with netCDF4.Dataset(filename, 'a') as nc:
            nc.delncattr('src_element_order')
            nc.delncattr('src_num_dofs_per_cell')
        rg3, sg3, dg3 = createRegridder()
        rg3.loadWeights(filename)
        assert numpy.array_equal(regrid(rg3, src_data), expected)

        # weights computed with another element are rejected
        with netCDF4.Dataset(filename, 'a') as nc:
            nc.src_element_order = numpy.int32(2)
        rg4, sg4, dg4 = createRegridder()
        with pytest.raises(RuntimeError):
            rg4.loadWeights(filename)
