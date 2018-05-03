from ugrid_reqder import UgridReader
from broken_line_iter import BrokenLineIter
from broken_segments_iter import BrokenSegmentsIter
import numpy


class RegridEdges:

	def __init__(self):
		
		# (dstCellEdge, srcCellEdge) -> weight
		self.weights = {}

		self.srcGrid = None
		self.srcLoc = None

		self.dstGrid = None


	def setSrcGridFile(self, filename):
		
		ur = UgridReader(filename)
		self.srcGrid = ur.getUnstructuredGrid()
		seld.srcLoc = ur.getUnstructuredGridLocator()


	def setDstGridFile(self, filename):
		
		ur = UgridReader(filename)
		self.dstGrid = ur.getUnstructuredGrid()


	def computeWeights(self):

        dstPtIds = vtk.vtkIdList()
		
        numSrcCells = self.srcGrid.GetNumberOfCells()
        numDstCells = self.dstGrid.GetNumberOfCells()

        # iterate over the dst grid cells
        for dstCellId in range(numDstCells):

            # iterate over the four edges of each dst cell
            self.dstGrid.GetCellPoints(dstCellId, dstPtIds)
            for i0 in range(4):
                i1 = (i0 + 1) % 4

                # get the start/end points of the dst edge
                dstVertId0 = dstPtIds.GetId(i0)
                dstVertId1 = dstPtIds.GetId(i1)
                dstEdgePt0 = self.dstGrid.GetPoint(dstVertId0)
                dstEdgePt1 = self.dstGrid.GetPoint(dstVertId1)

                # compute the intersections
                bli = BrokenLineIter([dstEdgePt0, dstEdgePt1])

                # find the intersection with the source grid
                bsi = BrokenSegmentsIter(self.srcGrid, self.srcLoc, bli)

                # compute the contribution to this edge
                for s in bsi:
                	srcCellId = s.getCellId()
                	k = (dstCellId, srcCellId)
                    xia = s.getBegCellParamCoord()
                    xib = s.getEndCellParamCoord()
                    self.weights[k] = self.weights.get(k, 0.) + basisIntegral(i0, xia, xib)


	def applyWeights(self, srcData):
		pass

