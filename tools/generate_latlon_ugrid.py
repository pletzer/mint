import netCDF4
import numpy
import sys
import argparse

"""
Generate grid and edge data on uniform grid and save result in UGRID file
"""

parser = argparse.ArgumentParser(description='Generate lat-lon grid and edge data in UGRID format')
parser.add_argument('-g', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-nx', default=1, type=int, 
                    help='Number of longitude cells')
parser.add_argument('-ny', default=1, type=int, 
                    help='Number of latitude cells')
parser.add_argument('-s', type=str, dest='stream_funct', default='sin(pi*x/180.)*cos(pi*y/180.)', help='Stream function of x (longitude in deg) and y (latitude in deg) used for setting the edge integrals')
args = parser.parse_args()

# check
if not args.grid_file:
    print("ERROR: must specify grid file (-g)")
    sys.exit(1)
try:
    grid_file, mesh_name = args.grid_file.split(':')
except:
    print("ERROR: grid file must be in the form 'FILENAME:MESHNAME'")
    sys.exit(2)


nx, ny = args.nx, args.ny

nc = netCDF4.Dataset(grid_file, 'w')

nnodes = (nx + 1) * (ny + 1)
nedges = nx * (ny + 1) + (nx + 1) * ny
nfaces = nx * ny
nnodesId = nc.createDimension('nnodes', nnodes)
nedgesId = nc.createDimension('nedges', nedges)
nfacesId = nc.createDimension('nfaces', nfaces)
fourId = nc.createDimension('four', 4)
twoId = nc.createDimension('two', 2)

mesh = nc.createVariable(mesh_name, "int", [])
mesh.cf_role = 'mesh_topology'
mesh.topology_dimension = 2
mesh.node_coordinates = 'lon lat'
mesh.face_node_connectivity = 'face_node'
mesh.face_edge_connectivity = 'face_edge'
mesh.edge_node_connectivity = 'edge_node'


faceNodeConn = nc.createVariable("face_node", "int64", ("nfaces", "four"))
faceNodeConn.cf_role = "face_node_connectivity"
faceNodeConn.start_index = 0

faceEdgeConn = nc.createVariable("face_edge", "int64", ("nfaces", "four"))
faceEdgeConn.cf_role = "face_edge_connectivity"
faceEdgeConn.start_index = 0

edgeNodeConn = nc.createVariable("edge_node", "int64", ("nedges", "two"))
edgeNodeConn.cf_role = "node_edge_connectivity"
edgeNodeConn.start_index = 0

xvar = nc.createVariable("lon", "float64", ("nnodes",))
xvar.standard_name = "longitude"
xvar.units = "degrees_east"

yvar = nc.createVariable("lat", "float64", ("nnodes",))
yvar.standard_name = "latitude"
yvar.units = "degrees_north"

edge_integrated_vel = nc.createVariable('edge_integrated_vel', 'float64', ('nedges',))
edge_integrated_vel.mesh = mesh_name
edge_integrated_vel.location = 'edge'

# set the lats/lons
lats = numpy.zeros((nnodes,), numpy.float64)
lons = numpy.zeros((nnodes,), numpy.float64)
for j in range(ny + 1):
    for i in range(nx + 1):
        index = i + (nx + 1)*j
        lons[index] = i * 360.0/float(nx)
        lats[index] = -90.0 + j * 180./float(ny)
xvar[:] = lons
yvar[:] = lats

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

# edge-node connectivity
en = numpy.zeros((nedges, 2), numpy.int64)
data = numpy.zeros((nedges,), numpy.float64)
# x edges
count = 0
for j in range(ny + 1):
    for i in range(nx):
        i00 = i + 0 + (nx + 1)*(j + 0)
        i10 = i + 1 + (nx + 1)*(j + 0)
        en[count, :] = i00, i10

        x, y = lons[i00], lats[i00]
        s00 = eval(args.stream_funct)
        x, y = lons[i10], lats[i10]
        s10 = eval(args.stream_funct)
        data[count] = s10 - s00

        count += 1
# y edges
for j in range(ny):
    for i in range(nx + 1):
        i00 = i + 0 + (nx + 1)*(j + 0)
        i01 = i + 0 + (nx + 1)*(j + 1)
        en[count, :] = i00, i01

        x, y = lons[i00], lats[i00]
        s00 = eval(args.stream_funct)
        x, y = lons[i01], lats[i01]
        s01 = eval(args.stream_funct)
        data[count] = s01 - s00

        count += 1

edgeNodeConn[...] = en
edge_integrated_vel[:] = data

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



