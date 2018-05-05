from ugrid_reader import UgridReader
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


    def setSrcGridFile(self, filename):
        """
        Set source grid
        @param file name containing UGRID description of the grid
        """ 
        ur = UgridReader(filename)
        self.srcGrid = ur.getUnstructuredGrid()
        self.srcLoc = ur.getUnstructuredGridCellLocator()
        self.srcLonLatPts = ur.getLonLatPoints()


    def setDstGridFile(self, filename):
        """
        Set destination grid
        @param file name containing UGRID description of the grid
        """    
        ur = UgridReader(filename)
        self.dstGrid = ur.getUnstructuredGrid()
        self.dstLonLatPts = ur.getLonLatPoints()


    def getNumSrcCells(self):
        return self.srcGrid.GetNumberOfCells()


    def getNumDstCells(self):
        return self.dstGrid.GetNumberOfCells()


    def getSrcLonLat(self):
        xy = self.srcLonLatPts.reshape((self.getNumSrcCells(), 4, 2))
        return xy[..., 0], xy[..., 1]


    def getDstLonLat(self):
        xy = self.dstLonLatPts.reshape((self.getNumDstCells(), 4, 2))
        return xy[..., 0], xy[..., 1]


    def computeWeights(self):
        """
        Compute the interpolation weights
        """
        raise NotImplementedError, 'ERROR: not implemented. You need to call the derived class method.'


    def applyWeights(self, srcData):
        """
        Apply the interpolation weigths to the field
        @param srcData line integrals on the source grid edges, dimensioned (ncells, 4)
        @return line integrals on the destination grid
        """
        numDstCells = self.dstGrid.GetNumberOfCells()
        res = numpy.zeros((numDstCells, 4), numpy.float64)
        for k in self.weights:
            dstCellId, srcCellId = k
            res[dstCellId, :] = self.weights[k].dot(srcData[srcCellId, :])
        return res


    def saveDstLoopData(self, dstData, filename):
        """
        Save destination loop integrals cell by cell
        @param dstData array of size numDstCells * 4 
        @param filename file name
        """
        numDstCells = self.getNumDstCells()
        self.dstLoopData = dstData[:, 0] + dstData[:, 1] - dstData[:, 2] - dstData[:, 3]
        self.dstData = vtk.vtkDoubleArray()
        self.dstData.SetNumberOfComponents(1)
        self.dstData.SetNumberOfTuples(numDstCells)
        self.dstData.SetVoidArray(self.dstLoopData, numDstCells, 1)
        self.dstGrid.GetCellData().SetScalars(self.dstData)

        writer= vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.dstGrid)
        writer.Update()




