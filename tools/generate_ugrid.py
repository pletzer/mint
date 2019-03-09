import netCDF4
import numpy


nx, ny = 1, 1

nc = netCDF4.Dataset('tiny.nc', 'w')

nnodes = (nx + 1) * (ny + 1)
nedges = nx * (ny + 1) + (nx + 1) * ny
nfaces = nx * ny
nnodesId = nc.createDimension('nnodes', nnodes)
nedgesId = nc.createDimension('nedges', nedges)
nfacesId = nc.createDimension('nfaces', nfaces)
fourId = nc.createDimension('four', 4)
twoId = nc.createDimension('two', 2)

faceNodeConn = nc.createVariable("face_node", "int64", ("nfaces", "four"))
faceNodeConn.cf_role = "face_node_connectivity"
faceNodeConn.start_index = 0
faceEdgeConn = nc.createVariable("face_edge", "int64", ("nfaces", "four"))
faceEdgeConn.cf_role = "face_edge_connectivity"
faceEdgeConn.start_index = 0
edgeNodeConn = nc.createVariable("edge_node", "int64", ("nedges", "two"))
edgeNodeConn.cf_role = "node_edge_connectivity"
edgeNodeConn.start_index = 0

x = nc.createVariable("lon", "float64", ("nnodes",))
x.standard_name = "longitude"
x.units = "degrees_east"
y = nc.createVariable("lat", "float64", ("nnodes",))
y.standard_name = "latitude"
y.units = "degrees_north"

# set the lats/lons
lats = numpy.zeros((nnodes,), numpy.float64)
lons = numpy.zeros((nnodes,), numpy.float64)
for j in range(ny + 1):
	for i in range(nx + 1):
		index = i + (nx + 1)*j
		lons[index] = i * 360.0/float(nx)
		lats[index] = -90.0 + j * 180./float(ny)

x[:] = lons
y[:] = lats


# face-node connectivity
fn = numpy.zeros((nfaces, 4), numpy.int64)
count = 0
for j in range(ny):
	for i in range(nx):
		i00 = i + 0 + (nx + 1)*(j + 0)
		i10 = i + 1 + (nx + 1)*(j + 0)
		i11 = i + 1 + (nx + 1)*(j + 1)
		i01 = i + 0 + (nx + 1)*(j + 1)
		fn[count, :] = i00, i10, i11, i01
		count += 1

faceNodeConn[...] = fn

en = numpy.zeros((nedges, 2), numpy.int64)

# edge-node connectivity
# x edges
count = 0
for j in range(ny + 1):
	for i in range(nx):
		i00 = i + 0 + (nx + 1)*(j + 0)
		i10 = i + 1 + (nx + 1)*(j + 0)
		en[count, :] = i00, i10
		count += 1

# y edges
for j in range(ny):
	for i in range(nx + 1):
		i00 = i + 0 + (nx + 1)*(j + 0)
		i01 = i + 0 + (nx + 1)*(j + 1)
		en[count, :] = i00, i01
		count += 1

edgeNodeConn[...] = en

# face-edge connectivity
fe = numpy.zeros((nfaces, 4), numpy.int64)
count = 0
for j in range(ny):
	for i in range(nx):
		is0 = 0 + i + 0 + (nx + 0)*(j + 0)
		is1 = 0 + i + 0 + (nx + 0)*(j + 1)
		i0s = nx*(ny + 1) + i + 0 + (nx + 1)*(j + 0)
		i1s = nx*(ny + 1) + i + 1 + (nx + 1)*(j + 0)
		fe[count, :] = is0, i1s, is1, i0s
		count += 1

faceEdgeConn[...] = fe


