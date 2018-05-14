import vtk
import numpy

class EdgeField:

	def __init__(self, grid):
		"""
		Constructor 
		@param grid instance of VTK unstructured grid
		"""
        self.pointArray = []
        self.pointData = vtk.vtkDoubleArray()
        self.points = vtk.vtkPoints()  
        self.grid = vtk.vtkUnstructuredGrid()
        self.loc = vtk.vtkCellLocator()


    def getLonLatPoints(self):
        """
        Get the longitudes and latitudes in radian at the cell vertices
        @return array
        """
        return self.pointArray[:, :2]

        
    def loadFromVtkFile(self, filename):
        """
        Load the grid and gfields form a VTK file
        @param filename VTK file
        """
        self.reader = vtk.vtkUnstructuredGridReader()
        self.reader.SetFileName(filename)
        self.reader.Update()
        self.grid = reader.GetOutput()
        self.locator.SetDataSet(self.grid)
        self.locator.BuildLocator()


    def saveToVtkFile(self, filename):
        """
        Save the grid to a VTK file
        @param filename VTK file
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.grid)
        writer.Update()



    def getUnstructuredGrid(self):
        """
        Get the unstructured grid
        @return vtkUnstructuredGrid instance
        """
        return self.grid


    def getUnstructuredGridCellLocator(self):
        """
        Get the unstructured grid cell locator
        @return vtkCellLocator instance
        """     
        return self.loc


    def setEdgeIntegralsFromStreamFunction(name, streamFuncData):
        """
        Compute edge integrals from a stream function given as a vertex field
        defined cell by cell 
        @param name name of the edge field
        @param streamFuncData arrays of size (numCells, 4)
        """
        # four edges per cell
        edgeVel = numpy.zeros(streamFuncData.shape, numpy.float64)
        for i0 in range(4):
            # edge direction is counter-clockwise
            i1 = (i0 + 1) % 4
            edgeVel[:, i0] = streamFuncData[:, i1] - streamFuncData[:, i0]

	def setFromStreamlineFunction(self, name, expr):
		"""
		Set field from expression
		@param name name of the field
        @param expr 
		"""
		psi = eval(expr)

	def getData(self):
		"""
		Get field as a data array
		"""
        pts = self.grid.GetPoints()
        x = pts.GetVoidPointer
		return self.dataArray
