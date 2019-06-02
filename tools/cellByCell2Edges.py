import vtk
import argparse

"""
Convert an edge field defined cell by cell to an edge field defined on a collection of line cells
"""

parser = argparse.ArgumentParser(description='Convert cell-by-cell edge field into an edge field')
parser.add_argument('-i', dest='inputFile', default='res.vtk', help='Specify path to VTK, cell-by-cell input file')
parser.add_argument('-o', dest='outputFile', default='edges.vtk', help='Specify name of VTK output file')
parser.add_argument('-v', dest='edgeFieldName', default='edge_integrated_velocity', help='Specify name of edge integrated variable')
args = parser.parse_args()

# read the file 
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.inputFile)
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
edata.SetName(args.edgeFieldName)
egrid.Allocate(nedges, 1)
data = ugrid.GetCellData().GetAbstractArray(args.edgeFieldName)
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
writer.SetFileName(args.outputFile)
writer.SetInputData(egrid)
writer.Update()

