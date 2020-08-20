import argparse
import sys
import numpy
import netCDF4

"""
Convert grid from spherical to cartesian
"""

LON_INDEX, LAT_INDEX = 0, 1
EPS = 1.234e-12

parser = argparse.ArgumentParser(description='Convert grid from spherical to cartesian')
parser.add_argument('-i', dest='inputFile', default='', help='Specify path to input Ugrid2D netcdf file (FILENAME:MESHANME)')
parser.add_argument('-o', dest='outputFile', default='', help='Specify path to output Ugrid2D netcdf file (FILENAME)')
args = parser.parse_args()

if len(args.inputFile) == 0:
    print('ERROR: must provide input Ugrid2D netcdf file (-i)')
    sys.exit(1)

if len(args.outputFile) == 0:
    print('ERROR: must provide output Ugrid2D netcdf file (-o)')
    sys.exit(2)

try:
    grid_file, grid_var = args.inputFile.split(':')
except:
    print('ERROR: could not extract grid name, specify -i FILENAME:MESHNAME')
    sys.exit(3)


def copyAttributes(fromNcVar, toNcVar):
    """
    Copy the attributes from fromNcVar to toNcVar
    """
    for attrName in fromNcVar.ncattrs():
        attrVal = getattr(fromNcVar, attrName)
        setattr(toNcVar, attrName, attrVal)


ncIn = netCDF4.Dataset(grid_file, 'r')
ncOut = netCDF4.Dataset(args.outputFile, 'w')


try:
    ncGridVar = ncIn[grid_var]
except:
    print('ERROR: could not extract grid name {} from file {}'.format(grid_var, grid_file))
    print('possible grid names are:')
    count = 0
    for v in ncIn.variables:
        var = ncIn[v]
        if hasattr(var, 'cf_role') and var.cf_role == 'mesh_topology':
            print(v)
            count += 1
    print('found {} variables with attribute "cf_role" == "mesh_topology"'.format(count))
    sys.exit(4)

# get the vertex coordinate names
xName, yName = ncIn[grid_var].node_coordinates.split(' ')
edgeNodeConnName = ncIn[grid_var].edge_node_connectivity

# get grid dimensions
numPoints = ncIn[xName].shape[0]
numEdges = ncIn[ncGridVar.edge_node_connectivity].shape[0]
numFaces = ncIn[ncGridVar.face_node_connectivity].shape[0]


# read the coordinates
lonIn = ncIn.variables[xName][:]
latIn = ncIn.variables[yName][:]

# read the edge-node connectivity
edgeNodeConnectivity = ncIn.variables[edgeNodeConnName][:]

# subtract base index
start_index = getattr(ncIn.variables[edgeNodeConnName], 'start_index', 0)
edgeNodeConnectivity -= start_index

numEdgesDimName = ncIn.variables[edgeNodeConnName].dimensions[0]


# create dimensions and variables
twoDimName = 'two'
fourDimName = 'four'
numPointsDimName = 'num_points'
numEdgesDimName = 'num_edges'
numFacesDimName = 'num_faces'
ncOut.createDimension(twoDimName, 2)
ncOut.createDimension(fourDimName, 4)
ncOut.createDimension(numPointsDimName, numPoints)
ncOut.createDimension(numEdgesDimName, numEdges)
ncOut.createDimension(numFacesDimName, numFaces)

# copy the topology over
gridOutName = grid_var
gridVarOut = ncOut.createVariable(gridOutName, 'i4')

# save Outstream grid
for attrName in ncGridVar.ncattrs():
    attrVal = getattr(ncGridVar, attrName)
    setattr(gridVarOut, attrName, attrVal)

# new coordinates
xNameOut = 'xcart'
yNameOut = 'ycart'
zNameOut = 'zcart'
gridVarOut.node_coordinates = f'{xNameOut} {yNameOut} {zNameOut}'

faceNodeOut = ncOut.createVariable(gridVarOut.face_node_connectivity, 'i4', 
                               (numFacesDimName, fourDimName))
copyAttributes(ncIn[ncGridVar.face_node_connectivity], faceNodeOut)
faceNodeOut[:] = ncIn[ncGridVar.face_node_connectivity][:]

faceEdgeOut = ncOut.createVariable(gridVarOut.face_edge_connectivity, 'i4',
                               (numFacesDimName, fourDimName))
copyAttributes(ncIn[ncGridVar.face_edge_connectivity], faceEdgeOut)
faceEdgeOut[:] = ncIn[ncGridVar.face_edge_connectivity][:]


edgeNodeOut = ncOut.createVariable(gridVarOut.edge_node_connectivity, 'i4',
                               (numEdgesDimName, twoDimName))
copyAttributes(ncIn[ncGridVar.edge_node_connectivity], edgeNodeOut)
edgeNodeOut[:] = ncIn[ncGridVar.edge_node_connectivity][:]

# save the new coordinates
xVarOut = ncOut.createVariable(xNameOut, 'f8', (numPointsDimName,))
xVarOut.standard_name = 'x-coordinate in Cartesian system'

yVarOut = ncOut.createVariable(yNameOut, 'f8', (numPointsDimName,))
yVarOut.standard_name = 'y-coordinate in Cartesian system'

zVarOut = ncOut.createVariable(zNameOut, 'f8', (numPointsDimName,))
zVarOut.standard_name = 'z-coordinate in Cartesian system'

# write
xVarOut[:] = numpy.cos(latIn*numpy.pi/180.)*numpy.cos(lonIn*numpy.pi/180.)
yVarOut[:] = numpy.cos(latIn*numpy.pi/180.)*numpy.sin(lonIn*numpy.pi/180.)
zVarOut[:] = numpy.sin(latIn*numpy.pi/180.)


ncIn.close()

ncOut.close()

