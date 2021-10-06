from mint import RegridEdges
import numpy
from pathlib import Path


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')


def test_compute_weights():

    # create a regridder
    rg = RegridEdges()

    # src is lon-lat grid so set the flags to 0, 0. This will
    # not apply any transformations to the grid
    rg.setSrcGridFlags(0, 0)

    # src grid file
    src_file = str(DATA_DIR / Path('latlon4x2.nc'))
    # load the src grid, latlon is the name of the mesh stored in
    # the netcdf file
    rg.loadSrcGrid(f'{src_file}$latlon')

    # dst is cubed-sphere. Cells at the poles need to be fixed to
    # make them as compact as possible. Use flags 1, 1
    rg.setDstGridFlags(1, 1)

    # dst grid file
    dst_file = str(DATA_DIR / Path('cs_4.nc'))
    # load the dst grid, physics is the name of the mesh stored in
    # the netcdf file
    rg.loadDstGrid(f'{dst_file}$physics')

    # compute the regridding weights. numCellsPerBucket is used to
    # accelerate the cell search. The smaller numCellPerBucket the
    # faster the search. However, there are edge cases where the
    # search fails when numCellsPerBucket is too small. periodX is
    # the periodicity length to add/subtract to make the cells well
    # behaved (periodX can be 0 if a regional model)
    rg.build(numCellsPerBucket=128, periodX=360., debug=2)

    # save the weights in a netCDF file
    rg.dumpWeights('test_regrid_edges_py.nc')


def test_apply_weights():
    # create a regridder
    rg = RegridEdges()

    # src is lon-lat grid
    rg.setSrcGridFlags(0, 0)
    rg.loadSrcGrid(f'{DATA_DIR}/latlon4x2.nc$latlon')

    # dst is cubed-sphere
    rg.setDstGridFlags(1, 1)
    rg.loadDstGrid(f'{DATA_DIR}/cs_4.nc$physics')

    # load the weights
    rg.loadWeights('test_regrid_edges_py.nc')

    num_src_edges = rg.getNumSrcEdges()
    num_dst_edges = rg.getNumDstEdges()
    print(f'number of edges (src/dst): {num_src_edges}/{num_dst_edges}')

    # create some mock field
    src_data = numpy.array(range(0, num_src_edges), numpy.float64)

    # allocate some space to receive the data
    dst_data = numpy.empty((num_dst_edges,), numpy.float64)

    # apply the weights to the src field, will fill in dst_data
    rg.apply(src_data, dst_data)

    check_sum = numpy.abs(dst_data).sum()
    assert(abs(check_sum - 515.8441563902018) < 1.e-8)
    print(f'check sum test passsed: {check_sum}')


if __name__ == '__main__':

    test_compute_weights()
    test_apply_weights()
