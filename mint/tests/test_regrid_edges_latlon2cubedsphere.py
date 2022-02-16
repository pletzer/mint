import mint
import numpy
import math
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

def streamFunc(point):

    x, y, _ = point
    return numpy.sin(x*numpy.pi/180.) * numpy.cos(y*numpy.pi/180.)


def test_1():

    # source grid
    sgrid = mint.Grid()
    # lat-lon
    sgrid.setFlags(0, 0)
    sgrid.loadFromUgrid2D(f'{DATA_DIR}/latlon100x50.nc$latlon')

    # destination grid
    dgrid = mint.Grid()
    # cubed-sphere
    dgrid.setFlags(1, 1)
    dgrid.loadFromUgrid2D(f'{DATA_DIR}/lfric_diag_wind.nc$Mesh2d')

    regridder = mint.RegridEdges()
    regridder.setSrcGrid(sgrid)
    regridder.setDstGrid(dgrid)
    regridder.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=False)
    regridder.computeWeights(debug=2)

    # get the cell-by-cell points
    spoints = sgrid.getPoints()

    # create src data
    numSrcCells = sgrid.getNumberOfCells()
    sdata = numpy.zeros((numSrcCells, 4), numpy.float64)

    # allocate dst data array
    numDstCells = dgrid.getNumberOfCells()
    ddata = numpy.zeros((numDstCells, 4), numpy.float64)

    for icell in range(numSrcCells):
        for i0 in range(4):
            # i1 is the second point of the edge in
            # anticlockwise direction
            i1 = (i0 + 1) % 4
            p0 = spoints[icell, i0, :]
            p1 = spoints[icell, i1, :]
            # our convention is to have the edges pointing in the
            # positive parametric direction, the last two edges
            # have the wrong sign
            sign = 1 - 2*(i0 // 2)
            # associate the same vorticity to each edge
            sdata[icell, i0] = sign * (streamFunc(p1) - streamFunc(p0))

    # apply the weights
    regridder.apply(sdata, ddata, placement=mint.CELL_BY_CELL_DATA)

    # compute the flow across a broken line
    xyz = numpy.array([(-170., -80., 0.), (170., 80., 0.)])

    sflux = mint.PolylineIntegral()
    sflux.setGrid(sgrid)
    sflux.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=False)
    sflux.computeWeights(xyz)

    dflux = mint.PolylineIntegral()
    dflux.setGrid(dgrid)
    dflux.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=False)
    dflux.computeWeights(xyz)

    exactFlux = streamFunc(xyz[-1, :]) - streamFunc(xyz[0, :])
    srcFluxVal = sflux.getIntegral(sdata)
    errSrcFlux = srcFluxVal - exactFlux
    dstFluxVal = dflux.getIntegral(ddata)
    errDstFlux = dstFluxVal - exactFlux
    print(f'fluxes: exact={exactFlux} over grid src: {srcFluxVal} (error={errSrcFlux:.2g}) dst: {dstFluxVal} (error={errDstFlux:.2g})')

    # save the grids and fields to VTK files
    sgrid.dump('sgrid.vtk')
    dgrid.dump('dgrid.vtk')

    mint.printLogMessages()
    mint.writeLogMessages('test_regrid_edges_latlon2cubedsphere.log')

if __name__ == '__main__':
    test_1()