import mint
import numpy

def main():

	# source grid
	sgrid = mint.Grid()
	# cubed-sphere
	sgrid.setFlags(1, 1)
	sgrid.loadFromUgrid2D('../../data/lfric_diag_wind.nc$Mesh2d')

	# destination grid
	dgrid = mint.Grid()
	# lat-lon
	dgrid.setFlags(0, 0)
	dgrid.loadFromUgrid2D('../../data/latlon100x50.nc$latlon')

	regridder = mint.RegridEdges()
	regridder.setSrcGrid(sgrid)
	regridder.setDstGrid(dgrid)
	regridder.buildLocator(numCellsPerBucket=128, periodX=360.0, enableFolding=0)
	regridder.computeWeights(debug=2)

	



if __name__ == '__main__':
	main()