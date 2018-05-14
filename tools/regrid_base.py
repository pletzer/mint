from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
import numpy
import vtk
import ctypes


class RegridBase(object):
    """
    Base class for cell by cell regridding
    """

    ZERO4x4 = numpy.zeros((4,4), numpy.float64)

    def __init__(self):
        """
        Constructor
        no args
        """
        
        # (dstCellId, srcCellId) -> weight as a 4x4 matrix
        self.weights = {}

        self.srcReader = vtk.vtkUnstructuredGridReader()
        self.srcGrid = vtk.vtkUnstructuredGrid()
        self.srcLoc = vtk.vtkCellLocator()
        self.srcEdgeField = []

        self.dstReader = vtk.vtkUnstructuredGridReader()
        self.dstGrid = vtk.vtkUnstructuredGrid()


    def setSrcFile(self, filename):
        """
        Set source grid and data from file
        @param reader instance of ReaderBase
        """
        self.srcReader.SetFileName(filename)
        self.srcReader.Update()
        self.srcGrid = self.srcReader.GetOutput()

        self.srcLoc.SetDataSet(self.srcGrid)
        self.srcLoc.BuildLocator()


    def setDstFile(self, filename):
        """
        Set destination grid from file
        @param reader instance of ReaderBase
        """    
        self.dstReader.SetFileName(filename)
        self.dstReader.Update()
        self.dstGrid = self.dstReader.GetOutput()


    def getSrcLonLat(self):
        """
        Get the source grid longitudes and latitudes in rad
        @return lon, lat arrays
        """
        return self.__getGridLonLat(self.srcGrid)


    def getDstLonLat(self):
        """
        Get the source grid longitudes and latitudes in rad
        @return lon, lat arrays
        """
        return self.__getGridLonLat(self.dstGrid)


    def computeWeights(self):
        """
        Compute the interpolation weights
        """
        raise NotImplementedError, 'ERROR: not implemented. You need to call the derived class method.'


    def applyWeights(self, srcData):
        """
        Apply the interpolation weights to the field
        @param srcData line integrals on the source grid edges, dimensioned (numSrcCells, 4)
        @return line integrals on the destination grid, array dimensioned (numDstCells, 4)
        """
        numDstCells = self.dstGrid.GetNumberOfCells()
        # allocate space
        res = numpy.zeros((numDstCells, 4), numpy.float64)
        # itererate over the desCellId, srcCellId pairs
        for k in self.weights:
            dstCellId, srcCellId = k
            res[dstCellId, :] += self.weights[k].dot(srcData[srcCellId, :])
        return res


    def getSrcEdgeData(self, name):
        """
        Get source edge data by name
        @param name name of the array
        @return numpy array
        """
        cellData = self.srcGrid.GetCellData().GetAbstractArray(name)
        numCells = self.srcGrid.GetNumberOfCells()
        return self.__getNumpyArrayFromVtkDoubleArray(numCells, 1, cellData)


    def saveDstLoopData(self, dstData, filename):
        """
        Save destination loop integrals cell by cell
        @param dstData array of size numDstCells * 4 
        @param filename file name
        """
        numDstCells = self.dstGrid.GetNumberOfCells()
        self.dstLoopData = dstData.sum(axis=1)
        self.dstData = vtk.vtkDoubleArray()
        self.dstData.SetNumberOfComponents(1)
        self.dstData.SetNumberOfTuples(numDstCells)
        self.dstData.SetVoidArray(self.dstLoopData, numDstCells, 1)
        self.dstGrid.GetCellData().SetScalars(self.dstData)

        writer= vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.dstGrid)
        writer.Update()


    def __getNumpyArrayFromVtkDoubleArray(self, numCells, numComponents, vtkArray):
        """
        Get a numpy array from a VTK array
        @param numCells number of cells
        @param numComponents number of components
        @param vtkArray instance of vtkDoubleArray
        @return numpy array
        """
        # vtkArray.GetVoidPointer(0) returns the address of the pointer as a string
        # "_X_void_p" where X is the hex address. Then convert X
        # to an int using base 16.
        addr = int(vtkArray.GetVoidPointer(0).split('_')[1], 16)
        ArrayType = ctypes.c_double * numCells * 4 * numComponents
        if numComponents == 1:
            return numpy.frombuffer(ArrayType.from_address(addr)).reshape((numCells, 4))
        else:
            return numpy.frombuffer(ArrayType.from_address(addr)).reshape((numCells, 4, numComponents))


    def __getGridLonLat(self, grid):
        """
        Get the grid longitudes and latitudes in rad
        @param grid instance of vtkUnstructuredGrid
        @return lon, lat arrays
        """
        data = grid.GetPoints().GetData()
        numCells = grid.GetNumberOfCells()
        xyz = self.__getNumpyArrayFromVtkDoubleArray(numCells, 3, data)
        lons, lats = xyz[..., 0], xyz[..., 1]
        return lons, lats


