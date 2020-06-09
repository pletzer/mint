import netCDF4
import numpy
import sys
import argparse

from numpy import sin, cos, pi, heaviside, exp

"""
Generate grid and edge data on uniform grid and save result in UGRID file
"""

parser = argparse.ArgumentParser(description='Generate output file storing the lat-lon grid and edge data in UGRID format')
parser.add_argument('-o', dest='file_mesh', default='', 
                    help='Specify the netcdf file containing the grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-nx', default=1, type=int, 
                    help='Number of longitude cells')
parser.add_argument('-ny', default=1, type=int, 
                    help='Number of latitude cells')
parser.add_argument('-nz', default=1, type=int, 
                    help='Number of elevation cells')
parser.add_argument('-nt', default=1, type=int, 
                    help='Number of time cells')
parser.add_argument('-s', type=str, dest='potential_funct', default='sin(pi*x/180.)*cos(pi*y/180.)*z**2 * exp(t)',
                    help='Potential function of x (lon in deg), y (lat in deg), z (0 <= z <=1) and t (t=0, ...nt-1)')
args = parser.parse_args()

# check
if not args.file_mesh:
    print("ERROR: must specify grid file (-o FILENAME:MESHNAME)")
    sys.exit(1)
try:
    grid_file, mesh_name = args.file_mesh.split(':')
except:
    print("ERROR: grid file must be in the form 'FILENAME:MESHNAME'")
    sys.exit(2)


nx, ny, nz, nt = args.nx, args.ny, args.nz, args.nt

nc = netCDF4.Dataset(grid_file, 'w')
nc.command = ' '.join(sys.argv)

nnodes = (nx + 1) * (ny + 1)
nedges = nx * (ny + 1) + (nx + 1) * ny
nfaces = nx * ny
nnodesId = nc.createDimension('nnodes', nnodes)
nedgesId = nc.createDimension('nedges', nedges)
nfacesId = nc.createDimension('nfaces', nfaces)
nelevsId = nc.createDimension('nelevs', nz + 1)
ntimesId = nc.createDimension('ntimes', nt + 1)
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

zvar = nc.createVariable("elev", "float64", ("nelevs",))
zvar.setncattr('name', "full_levels")

tvar = nc.createVariable("time", "float64", ("ntimes",))
tvar.standard_name = "time"
tvar.calendar = "gregorian"
tvar.units = "hours since 2016-01-01 15:00:00"
tvar.time_origin = "2016-01-01 15:00:00"

u = nc.createVariable('u', 'float64', ('ntimes', 'nelevs', 'nedges',))
u.mesh = mesh_name
u.location = 'edge'

potentialvar = nc.createVariable('potential_funct', 'float64', ('ntimes', 'nelevs', 'nnodes',))
potentialvar.mesh = mesh_name
potentialvar.location = 'node'

# set the lats/lons
lats = numpy.zeros((nnodes,), numpy.float64)
lons = numpy.zeros((nnodes,), numpy.float64)
dlat, dlon = 180./float(ny), 360.0/float(nx)
index = 0
for j in range(ny + 1):
    y = -90.0 + j*dlat
    for i in range(nx + 1):
        x = -180. + i*dlon
        lons[index] = x
        lats[index] = y
        index += 1
xvar[:] = lons
yvar[:] = lats

# set the potential values
point_data = numpy.zeros((nt + 1, nz + 1, nnodes,), numpy.float64)
for l in range(nt + 1):
    t = float(l)
    for k in range(nz + 1):
        # elevation
        z = float(k)
        index = 0
        for j in range(ny + 1):
            for i in range(nx + 1):
                x, y = lons[index], lats[index]
                point_data[l, k, index] = eval(args.potential_funct)
                index += 1
potentialvar[:] = point_data

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
edge_data = numpy.zeros((nt + 1, nz + 1, nedges,), numpy.float64)

# x edges
for l in range(nt + 1):
    t = float(l)
    for k in range(nz + 1):
        # elevation
        z = float(k)

        count = 0

        # x edges
        for j in range(ny + 1):
            for i in range(nx):
                i00 = i + 0 + (nx + 1)*(j + 0)
                i10 = i + 1 + (nx + 1)*(j + 0)
                en[count, :] = i00, i10

                x, y = lons[i00], lats[i00]
                s00 = eval(args.potential_funct)
                x, y = lons[i10], lats[i10]
                s10 = eval(args.potential_funct)
                edge_data[l, k, count] = s10 - s00

                count += 1

        # y edges
        for j in range(ny):
            for i in range(nx + 1):
                i00 = i + 0 + (nx + 1)*(j + 0)
                i01 = i + 0 + (nx + 1)*(j + 1)
                en[count, :] = i00, i01

                x, y = lons[i00], lats[i00]
                s00 = eval(args.potential_funct)
                x, y = lons[i01], lats[i01]
                s01 = eval(args.potential_funct)
                edge_data[l, k, count] = s01 - s00

                count += 1

edgeNodeConn[...] = en
u[:] = edge_data

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


nc.close()



