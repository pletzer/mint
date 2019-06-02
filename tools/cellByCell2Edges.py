import vtk

inputFile = '../build2/res.vtk'
outputFile = 'edges.vtk'
edgeFieldName = 'edge_integrated_velocity'

# read the file 
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(inputFile)
reader.Update()

# get the unstructured grid
ugrid = reader.GetOutput()

# create a new unstructured grid made of edges only
edata = vtk.vtkDoubleArray()
egrid = vtk.vtkUnstructuredGrid()
egrid.SetPoints(ugrid.GetPoints())
ncells = ugrid.GetNumberOfCells()
nedges = 4 * ncells
edata.SetNumberOfComponents(1)
edata.SetNumberOfTuples(nedges)
edata.SetName(edgeFieldName)
egrid.Allocate(nedges, 1)
data = ugrid.GetCellData().GetAbstractArray(edgeFieldName)
iedge = 0
for icell in range(ncells):
	cell = ugrid.GetCell(icell)
	nedgesPerCell = cell.GetNumberOfEdges()
	for ie in range(nedgesPerCell):
		edge = cell.GetEdge(ie)
		edgePtIds = edge.GetPointIds()
		egrid.InsertNextCell(vtk.VTK_LINE, edgePtIds)
		edata.SetTuple(iedge, (data.GetComponent(icell, ie),))
		iedge += 1

# add the adat to the edge grid
egrid.GetCellData().AddArray(edata)

# write the edge grid
writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName(outputFile)
writer.SetInputData(egrid)
writer.Update()

