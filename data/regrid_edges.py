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
        # iterate over the dst grid edges
        for idstC in range(numDstCells):
            # iterate over the four edges of each dst cell
            self.dstGrid.GetCellPoints(idstC, dstPtIds)
            for i0 in range(4):
                i1 = (i0 + 1) % 4
                # get the start/end points of the dst edge
                dstVertId0 = dstPtIds.GetId(i0)
                dstVertId1 = dstPtIds.GetId(i1)
                p0 = self.dstGrid.GetPoint(dstVertId0)
                p1 = self.dstGrid.GetPoint(dstVertId1)

                # compute the intersections
                bli = BrokenLineIter([p0, p1])
                bsi = BrokenSegmentsIter(self.srcGrid, self.srcLoc, bli)

                # compute the contribution to this edge
                for s in bsi:
                    pass



	def applyWeights(self, srcData):
		pass

