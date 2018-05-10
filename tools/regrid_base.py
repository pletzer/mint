from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
import numpy
import vtk

    
def edgeIntegralFromStreamFunction(streamFuncData):
    """
    Compute edge integrals from a stream function given as a vertex field
    defined cell by cell 
    @param streamFuncData arrays of size (numCells, 4)
    @return edge data
    """
    # four edges per cell
    edgeVel = numpy.zeros(streamFuncData.shape, numpy.float64)
    for i0 in range(4):
        # edge direction is counter-clockwise
        i1 = (i0 + 1) % 4
        edgeVel[:, i0] = streamFuncData[:, i1] - streamFuncData[:, i0]
    return edgeVel



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


    def setSrcGridFile(self, filename, format='ugrid', padding=0):
        """
        Set source grid
        @param file name containing UGRID description of the grid
        @param format either 'ugrid' or 'um'
        @param padding number of additional lon cells to add to the high end
                       (only for grid in UM format)
        """ 
        ur = None
        if format == 'ugrid':
            ur = UgridReader(filename)
        else:
            ur = LatLonReader(filename, padding=padding)
        self.srcGrid = ur.getUnstructuredGrid()
        self.srcLoc = ur.getUnstructuredGridCellLocator()
        self.srcLonLatPts = ur.getLonLatPoints()


    def setDstGridFile(self, filename, format='ugrid'):
        """
        Set destination grid
        @param file name containing UGRID description of the grid
        """    
        ur = None
        if format == 'ugrid':
            ur = UgridReader(filename)
        else:
            ur = LatLonReader(filename)
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
        numDstCells = self.getNumDstCells()
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



