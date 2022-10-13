import netCDF4
import defopt
import numpy

def main(*, nx: int=10, ny: int=5, output_file: str='sgrid.nc'):
	"""
	Generate sgrid vector data
	:param nx: number of x cells
	:param ny: number of y cells

	"""

	# vertex grid
	nx1 = nx + 1
	ny1 = ny + 1
	x1d = numpy.linspace(-180., 180, nx1)
	y1d = numpy.linspace(-90., 90., ny1)
	xx, yy = numpy.meshgrid(x1d, y1d, indexing='xy')

	# u points
	xxu = xx
	yyu = 0.5*(yy[:-1, :] + yy[1:, :])

	# v points
	xxv = 0.5*(xx[:, :-1] + xx[:, 1:])
	yyv = yy

	deg2rad = numpy.pi/180.
	uu = numpy.cos(deg2rad * yyu)
	vv = numpy.sin(deg2rad * xxv)

	# write the data to file
	with netCDF4.Dataset(output_file, 'w') as nc:

		# dimensions
		nx1dim = nc.createDimension('nx1', nx1)
		ny1dim = nc.createDimension('ny1', ny1)
		nxdim = nc.createDimension('nx', nx)
		nydim = nc.createDimension('ny', ny)
		nxudim = nc.createDimension('nxu', uu.shape[1])
		nyudim = nc.createDimension('nyu', uu.shape[0])
		nxvdim = nc.createDimension('nxv', vv.shape[1])
		nyvdim = nc.createDimension('nyv', vv.shape[0])

		xvar = nc.createVariable('x', 'f4', ('ny1', 'nx1'))
		xvar.standard_name = 'longitude'
		xvar.units = 'degree_east'

		yvar = nc.createVariable('y', 'f4', ('ny1', 'nx1'))
		yvar.standard_name = 'latitude'
		yvar.units = 'degree_north'

		uvar = nc.createVariable('u', 'f4', ('nyu', 'nxu'))
		uvar.description = 'x-velocity'
		uvar.units = 'm s-1'
		uvar.grid = 'grid'
		uvar.location = 'edge1' # Arakawa C

		vvar = nc.createVariable('v', 'f4', ('nyv', 'nxv'))
		vvar.description = 'y-velocity'
		vvar.units = 'm s-1'
		vvar.grid = 'grid'
		vvar.location = 'edge2' # Arakawa C

		gridvar = nc.createVariable('grid', 'int32', [])
		gridvar.cf_role = 'grid_topology'
		gridvar.topology_dimensions = numpy.int32(2)
		gridvar.node_dimensions = "nx1 ny1"
		gridvar.node_coordinates = "x y"
		gridvar.face_dimensions = "nx: nx1 (padding: node) ny: ny1 (padding: node)"

		# write the data
		xvar[:] = xx
		yvar[:] = yy
		uvar[:] = uu
		vvar[:] = vv



if __name__ == '__main__':
	defopt.run(main)



