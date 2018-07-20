import netCDF4
import numpy
import vtk

class ReaderBase(object):

    PERIODICITY_LENGTH = 360.

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


    def saveToVtkFile(self, filename, binary=True):
        """
        Save the grid to a VTK file
        @param filename VTK file
        """
        writer = vtk.vtkUnstructuredGridWriter()
        if binary:
            writer.SetFileTypeToBinary()
        writer.SetFileName(filename)
        writer.SetInputData(self.vtk['grid'])
        writer.Update()



    def getUnstructuredGrid(self):
        """
        Get the unstructured grid
        @return vtkUnstructuredGrid instance
        """
        return self.vtk['grid']


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

    def setPointField(self, name, data):
        """
        Set point field
        @param name name of the field
        @param data array of size (numCells, 4)
        """
        self.pointArray = data
        self.pointData = vtk.vtkDoubleArray()
        self.pointData.SetName(name)
        self.pointData.SetNumberOfComponents(1)
        numCells = self.getNumberOfCells()
        self.pointData.SetNumberOfTuples(numCells*4)
        self.pointData.SetVoidArray(self.pointArray, numCells*4, 1)
        self.vtk['grid'].GetPointData().AddArray(self.pointData)        


    def setPointVectorField(self, name, uData, vData):
        """
        Set vector field on cell points
        @param name name of the field
        @param data array of size (numCells, 4)
        """
        numCells = self.getNumberOfCells()
        self.pointVectorArray = numpy.zeros((numCells*4, 3), numpy.float64)
        self.pointVectorArray[:, 0] = uData.flat
        self.pointVectorArray[:, 1] = vData.flat
        self.pointVectorData = vtk.vtkDoubleArray()
        self.pointVectorData.SetName(name)
        # 4 points per cell, 3 components
        self.pointVectorData.SetNumberOfComponents(3)
        self.pointVectorData.SetNumberOfTuples(numCells*4)
        self.pointVectorData.SetVoidArray(self.pointVectorArray, numCells*4*3, 1)
        self.vtk['grid'].GetPointData().AddArray(self.pointVectorData)        

    def getEdgeFieldFromStreamData(self, streamFuncData):
        """
        Get the edge integrated values from the nodal stream function data
        @param streamFuncData stream function data
        """
        numCells = self.getNumberOfCells()
        edgeArray = numpy.zeros((numCells, 4), numpy.float64)
        for i0 in range(4):
            # edge direction is counter-clockwise
            i1 = (i0 + 1) % 4
            edgeArray[:, i0] = streamFuncData[:, i1] - streamFuncData[:, i0]
        return edgeArray

    def getLoopIntegralsFromStreamData(self, streamFuncData):
        """
        Get the cell loop integral from nodal the stream function data
        @param streamFuncData stream function data
        """
        numCells = self.getNumberOfCells()
        cellArray = numpy.zeros((numCells,), numpy.float64)
        for i0 in range(4):
            # edge direction is counter-clockwise
            i1 = (i0 + 1) % 4
            cellArray[:] += streamFuncData[:, i1] - streamFuncData[:, i0]
        return cellArray

    def setLoopIntegrals(self, name, data):
        """
        Set cell loop integral field
        @param name name of the field
        @param data array of size (numCells,)
        """
        self.cellArray = data
        self.cellData = vtk.vtkDoubleArray()
        self.cellData.SetName(name)
        self.cellData.SetNumberOfComponents(1)
        numCells = self.getNumberOfCells()
        self.cellData.SetNumberOfTuples(numCells)
        self.cellData.SetVoidArray(self.cellArray, numCells*1, 1)
        self.vtk['grid'].GetCellData().AddArray(self.cellData)


