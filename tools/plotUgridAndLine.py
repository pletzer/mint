import vtk
import netCDF4

srcFile = '../data/latlon8x4.nc'
meshName = 'latlon'
targetPoints = [(0., 22.5, 0.), (22.5, 20.9410204722438422, 0.)]

nc = netCDF4.Dataset(srcFile)
mesh = nc.variables[meshName]
f2n = nc.variables[mesh.face_node_connectivity][:]
f2e = nc.variables[mesh.face_edge_connectivity][:]
e2n = nc.variables[mesh.edge_node_connectivity][:]

lons = nc.variables[mesh.node_coordinates.split()[0]][:]
lats = nc.variables[mesh.node_coordinates.split()[1]][:]

numPoints = len(lons)
numEdges = e2n.shape[0]
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
writer.SetFileName('latlon8x4.vtk')
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



