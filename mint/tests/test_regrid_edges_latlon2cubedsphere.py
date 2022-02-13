import mint
import numpy
from pathlib import Path

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

VORTEX_CENTRE = numpy.array([10., 20., 0.])
VORTEX_SIG = 30.0

def vortexFunc(point):
	dx = point - VORTEX_CENTRE
	# gaussian bump
	return numpy.exp( -dx.dot(dx)/(2.*VORTEX_SIG*VORTEX_SIG))


def test_1():

	# source grid
	sgrid = mint.Grid()
	# cubed-sphere
	sgrid.setFlags(1, 1)
	sgrid.loadFromUgrid2D(f'{DATA_DIR}/lfric_diag_wind.nc$Mesh2d')

	# destination grid
	dgrid = mint.Grid()
	# lat-lon
	dgrid.setFlags(0, 0)
	dgrid.loadFromUgrid2D(f'{DATA_DIR}/latlon100x50.nc$latlon')

	regridder = mint.RegridEdges()
	regridder.setSrcGrid(sgrid)
	regridder.setDstGrid(dgrid)
	regridder.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=0)
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
		# get the mid point of the cell
		midPoint = spoints[icell, :, :].mean(axis=0)
		# get the vortex strength for this cell
		vortexStrength = vortexFunc(midPoint)
		# iterate over the vertices of the cell,
		# i0 is the first point of the edge
		for i0 in range(4):
			# i1 is the second point of the edge in
			# anticlockwise direction
			i1 = (i0 + 1) % 4
			# our convention is to have the edges pointing in the
			# positive parametric direction, the last two edges
			# have the wrong sign
			sign = 1 - 2*(i0 // 2)
			# associate the same vorticity to each edge
			sdata[icell, i0] = sign * 0.25*vortexStrength

	# apply the weights
	regridder.apply(sdata, ddata, placement=mint.CELL_BY_CELL_DATA)

	# attach fields to the src and dst grids
	sgrid.attach('vorticity', sdata[:, 0] + sdata[:, 1] - sdata[:, 2] - sdata[:, 3])
	dgrid.attach('vorticity', ddata[:, 0] + ddata[:, 1] - ddata[:, 2] - ddata[:, 3])

	# save the grids and fields to VTK files
	sgrid.dump('sgrid.vtk')
	dgrid.dump('dgrid.vtk')

	mint.writeLogMessages('test_regrid_edges_latlon2cubedsphere.log')

if __name__ == '__main__':
	test_1()