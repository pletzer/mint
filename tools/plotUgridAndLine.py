import vtk
import netCDF4
import argparse
import re
import sys
import numpy

parser = argparse.ArgumentParser(description='Plot grid and line')
parser.add_argument('-i', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-p', dest='points', default='(0.0,0.0),(360.,0.0)', 
                    help='Specify the points of the line')
parser.add_argument('-e', dest='edge_field_name', default='edge_integrated_velocity', 
                    help='Specify the name of the edge integrated field')


args = parser.parse_args()


srcFile, meshName = args.grid_file.split(':')

# parse the target points
targetPoints = re.sub('\s+', '', args.points, )
targetPoints = re.sub(r'^\(', '', targetPoints)
targetPoints = re.sub(r'\)$', '', targetPoints)
targetPoints = [ eval(p + ', 0.0') for p in targetPoints.split('),(') ]

#[(0., 22.5, 0.), (22.5, 20.9410204722438422, 0.)]

# read the data
nc = netCDF4.Dataset(srcFile)
mesh = nc.variables[meshName]

f2n = nc.variables[mesh.face_node_connectivity][:]
f2n -= nc.variables[mesh.face_node_connectivity].start_index

f2e = nc.variables[mesh.face_edge_connectivity][:]
f2e -= nc.variables[mesh.face_edge_connectivity].start_index

e2n = nc.variables[mesh.edge_node_connectivity][:]
e2n -= nc.variables[mesh.edge_node_connectivity].start_index

lons = nc.variables[mesh.node_coordinates.split()[0]][:]
lats = nc.variables[mesh.node_coordinates.split()[1]][:]

edgeFieldVar = nc.variables[args.edge_field_name]
if not hasattr(edgeFieldVar, 'mesh'):
	print('ERROR: edge field ' + args.edge_field_name + ' must have attribute mesh')
	sys.exit(1)
if edgeFieldVar.mesh != meshName:
	print('ERROR: edge field ' + args.edge_field_name + "'s attribute mesh is not " + meshName)
	sys.exit(2)
if not hasattr(edgeFieldVar, 'location'):
	print('ERROR: edge field ' + args.edge_field_name + ' must have attribute location')
	sys.exit(3)
if edgeFieldVar.location != 'edge':
	print('ERROR: edge field ' + args.edge_field_name + "'s attribute location is not edge")
	sys.exit(4)
# read the edge field
edgeField = edgeFieldVar[:]	

numPoints = len(lons)
numFaces = f2n.shape[0]

# create the face grid
pointArray = vtk.vtkDoubleArray()
pointArray.SetNumberOfComponents(3)
pointArray.SetNumberOfTuples(numPoints)
for i in range(numPoints):
	pointArray.SetTuple(i, (lons[i], lats[i], 0.0))

points = vtk.vtkPoints()
points.SetData(pointArray)

gridFaces = vtk.vtkUnstructuredGrid()
gridFaces.SetPoints(points)
gridFaces.Allocate(1, 1)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(4)
edgeData = vtk.vtkDoubleArray()
edgeData.SetName(args.edge_field_name)
edgeData.SetNumberOfComponents(4)
edgeData.SetNumberOfTuples(numFaces)
for i in range(numFaces):

	# face to node connectivity, assume to be in counterclockwise direction
	counterClockwiseNodeIds = list(f2n[i, :])

	# set the cell points
	for j in range(4):
		ptIds.SetId(j, counterClockwiseNodeIds[j])

	gridFaces.InsertNextCell(vtk.VTK_QUAD, ptIds)

	# set the edge field values, edgeTuple is in the counterclockwise direction
	edgeTuple = [0., 0., 0., 0.]

	# get the edges attached to this face
	edgeIds = f2e[i, :]

	index = 0
	for edgeId in edgeIds:

		# get the nodes attached to this edge
		i0, i1 = e2n[edgeId, :]

		# get the index of the start node
		pos0 = counterClockwiseNodeIds.index(i0)

		# get the index of the end node
		pos1 = counterClockwiseNodeIds.index(i1)

		pLo = min(pos0, pos1)
		pHi = max(pos0, pos1)
		sign = 1
		if pLo // 2 == 0:
			# first two edges
			if pos1 < pos0:
				# clockwise direction
				sign = -1
		else:
			# last two edges
			if pos1 > pos0:
				# anticlockwise
				sign = -1

		edgeTuple[index] = sign * edgeField[edgeId]
		#print('face {} edge {} nodes {} {} pos0 {} pos1 {} sign {} pLo {} pHi {} anticlockwise nodes {}'.format(i, edgeId, i0, i1, pos0, pos1, sign, pLo, pHi, counterClockwiseNodeIds))
		index += 1

	print('loop integral for face {} is {}'.format(i, numpy.sum(edgeTuple[:2]) - numpy.sum(edgeTuple[2:])))
	edgeData.SetTuple(i, edgeTuple)


# add the edge field
gridFaces.GetCellData().AddArray(edgeData)

# write the grid file
writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(re.sub('\.nc', '.vtk', srcFile))
writer.SetInputData(gridFaces)
writer.Update()

# polyline grid
numPolyPoints = len(targetPoints)
polyPointArray = vtk.vtkDoubleArray()
polyPointArray.SetNumberOfComponents(3)
polyPointArray.SetNumberOfTuples(numPolyPoints)
for i in range(numPolyPoints):
	polyPointArray.SetTuple(i, targetPoints[i])

polyPoints = vtk.vtkPoints()
polyPoints.SetData(polyPointArray)

polyGrid = vtk.vtkUnstructuredGrid()
polyGrid.SetPoints(polyPoints)
polyGrid.Allocate(1, 1)
ptIds = vtk.vtkIdList()
ptIds.SetNumberOfIds(2)
for i in range(numPolyPoints - 1):
	ptIds.SetId(0, f2n[i, 0])
	ptIds.SetId(1, f2n[i, 1])
	polyGrid.InsertNextCell(vtk.VTK_LINE, ptIds)


# write the VTK file
writer.SetFileName('poly.vtk')
writer.SetInputData(polyGrid)
writer.Update()



