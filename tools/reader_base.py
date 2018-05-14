import netCDF4
import numpy
import vtk

class ReaderBase(object):

    TWOPI = 2. * numpy.pi

    def __init__(self):
        """
        Constructor
        No args
        """

        self.vtk = {
            'pointArray': [],
            'pointData': vtk.vtkDoubleArray(),
            'points': vtk.vtkPoints(),
            'grid': vtk.vtkUnstructuredGrid(),
            'locator': vtk.vtkCellLocator(),
        }


    def getLonLatPoints(self):
        """
        Get the longitudes and latitudes in radian at the cell vertices
        @return array
        """
        return self.vtk['pointArray'][:, :2]

        
    def loadFromVtkFile(self, filename):
        """
        Load the grid and gfields form a VTK file
        @param filename VTK file
        """
        self.reader = vtk.vtkUnstructuredGridReader()
        self.reader.SetFileName(filename)
        self.reader.Update()
        self.vtk['grid'] = reader.GetOutput()
        self.vtk['locator'].SetDataSet(self.vtk['grid'])
        self.vtk['locator'].BuildLocator()


    def saveToVtkFile(self, filename):
        """
        Save the grid to a VTK file
        @param filename VTK file
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.vtk['grid'])
        writer.Update()



    def getUnstructuredGrid(self):
        """
        Get the unstructured grid
        @return vtkUnstructuredGrid instance
        """
        return self.vtk['grid']


    def getUnstructuredGridCellLocator(self):
        """
        Get the unstructured grid cell locator
        @return vtkCellLocator instance
        """    	
    	return self.vtk['locator']

    def getNumberOfCells(self):
        """
        Get the number of cells 
        @return number
        """
        return self.vtk['grid'].GetNumberOfCells()


    def getLonLat(self):
        """
        Get the longitudes and latitudes as separate arrays
        @return lon and lat arrays of size (numCells, 4)
        """
        xy = self.vtk['pointArray'].reshape((self.getNumberOfCells(), 4, 3))
        return xy[..., 0], xy[..., 1]


    def setEdgeField(self, name, data):
        """
        Set edge field
        @param name name of the field
        @param data array of size (numCells, 4)
        """
        self.edgeArray = data
        self.edgeData = vtk.vtkDoubleArray()
        self.edgeData.SetName(name)
        self.edgeData.SetNumberOfComponents(4)
        numCells = self.getNumberOfCells()
        self.edgeData.SetNumberOfTuples(numCells)
        self.edgeData.SetVoidArray(self.edgeArray, numCells*4, 1)
        self.vtk['grid'].GetCellData().AddArray(self.edgeData)


    def getEdgeFieldFromStreamData(self, streamFuncData):
        """
        Set the edge integrated values from nodal stream function data
        @param streamFuncData stream function data
        """
        numCells = self.getNumberOfCells()
        edgeArray = numpy.zeros((numCells, 4), numpy.float64)
        for i0 in range(4):
            # edge direction is counter-clockwise
            i1 = (i0 + 1) % 4
            edgeArray[:, i0] = streamFuncData[:, i1] - streamFuncData[:, i0]
        return edgeArray


