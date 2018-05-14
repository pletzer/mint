from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
import numpy
import vtk


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

        self.srcGrid = None
        self.srcLoc = None
        self.srcLonLatPts = None

        self.dstGrid = None
        self.dstLonLatPts = None


    def setSrcGrid(self, reader):
        """
        Set source grid
        @param reader instance of ReraderBase
        """ 
        self.srcGrid = reader.getUnstructuredGrid()
        self.srcLoc = reader.getUnstructuredGridCellLocator()
        self.srcLonLatPts = reader.getLonLatPoints()


    def setDstGrid(self, reader):
        """
        Set destination grid
        @param reader instance of ReaderBase
        """    
        self.dstGrid = reader.getUnstructuredGrid()
        self.dstLonLatPts = reader.getLonLatPoints()


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
        # allocate space
        res = numpy.zeros((self.getNumDstCells(), 4), numpy.float64)
        # itererate over the desCellId, srcCellId pairs
        for k in self.weights:
            dstCellId, srcCellId = k
            res[dstCellId, :] += self.weights[k].dot(srcData[srcCellId, :])
        return res


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



