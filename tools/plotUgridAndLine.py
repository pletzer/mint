import vtk
import netCDF4
import argparse
import re

parser = argparse.ArgumentParser(description='Plot grid and line')
parser.add_argument('-i', dest='grid_file', default='', 
                    help='Specify the netcdf file containing the grid and the mesh name in the format "FILENAME:MESHNAME"')
parser.add_argument('-p', dest='points', default='(0.0,0.0),(360.,0.0)', 
                    help='Specify the points of the line')


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

lons = nc.variables[mesh.node_coordinates.split()[0]][:]
lats = nc.variables[mesh.node_coordinates.split()[1]][:]

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
for i in range(numFaces):
	ptIds.SetId(0, f2n[i, 0])
	ptIds.SetId(1, f2n[i, 1])
	ptIds.SetId(2, f2n[i, 2])
	ptIds.SetId(3, f2n[i, 3])
	gridFaces.InsertNextCell(vtk.VTK_QUAD, ptIds)

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

writer.SetFileName('poly.vtk')
writer.SetInputData(polyGrid)
writer.Update()



