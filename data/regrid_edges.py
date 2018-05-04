from regrid_base import RegridBase
from ugrid_reader import UgridReader
from broken_line_iter import BrokenLineIter
from broken_segments_iter import BrokenSegmentsIter
import numpy
import vtk
import ctypes


class RegridEdges(RegridBase):
    """
    Class for regridding edge field to edge field using a cell by cell approach
    """

    def __init__(self):
        """
        Constructor
        no args
        """
        super(RegridEdges, self).__init__()


    def computeWeights(self):
        """
        Compute the interpolation weights
        """

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
                for seg in bsi:
                    srcCellId = seg.getCellId()
                    xia = seg.getBegCellParamCoord()
                    xib = seg.getEndCellParamCoord()

                    dxi = xib - xia
                    xiMid = 0.5*(xia + xib)

                    k = (dstCellId, srcCellId)

                    # compute the weights from each src edge
                    ws = numpy.array([dxi[0]*(1.0 - xiMid[1]),
                                      dxi[1]*(0.0 + xiMid[0]),
                                      dxi[0]*(0.0 + xiMid[1]),
                                      dxi[1]*(1.0 - xiMid[0])])

                    if not self.weights.has_key(k):
                        # initialize the weights
                        self.weights[k] = self.ZERO4x4
                    self.weights[k][i0, :] += ws

                totalT = seg.getIntegratedParamCoord()
                if abs(totalT - 1.0) > 1.e-6:
                    print('Warning: total t of segment: {:.3f} != 1 (diff={:.1g}), dst cell {} point {} -> {}'.format(\
                        totalT, totalT - 1.0, dstCellId, dstEdgePt0, dstEdgePt1))


###############################################################################
def edgeIntegralFromStreamFunction(streamFuncData):
    edgeVel = numpy.zeros(streamFuncData.shape, numpy.float64)
    for i0 in range(4):
        i1 = (i0 + 1) % 4
        pm = 1 - 2*(i0 // 2) # + for i0 = 0, 1, - for i0 = 2, 3
        edgeVel[:, i0] = pm * (streamFuncData[:, i1] - streamFuncData[:, i0])
    return edgeVel


def main():
    from math import pi, sin, cos, log, exp
    import argparse

    parser = argparse.ArgumentParser(description='Regriod edge field')
    parser.add_argument('-s', dest='src', default='mesh_C4.nc', help='Specify UGRID source grid file')
    parser.add_argument('-d', dest='dst', default='mesh_C4.nc', help='Specify UGRID destination grid file')
    parser.add_argument('-S', dest='srcStreamFunc', default='x', 
                        help='Stream function as a function of x (longitude in rad) and y (latitude in rad)')
    args = parser.parse_args()

    rgrd = RegridEdges()
    rgrd.setSrcGridFile(args.src)
    rgrd.setDstGridFile(args.dst)
    rgrd.computeWeights()

    # compute stream function on cell vertices
    x, y = rgrd.getSrcLonLat()
    srcPsi = eval(args.srcStreamFunc)

    # compute edge integrals
    srcEdgeVel = edgeIntegralFromStreamFunction(srcPsi)

    # apply the weights 
    dstEdgeVel = rgrd.applyWeights(srcEdgeVel)

    # compute the exact edge field on the destination grid
    x, y = rgrd.getDstLonLat()
    dstPsi = eval(args.srcStreamFunc)
    dstEdgeVelExact = edgeIntegralFromStreamFunction(dstPsi)

    # compute the error
    diff = numpy.fabs(dstEdgeVelExact - dstEdgeVel)
    maxError = diff.max()
    minError = diff.min()
    print('Min/max error              : {}/{}'.format(minError, maxError))
    error = numpy.fabs(diff).sum()
    print('Sum of interpolation errors: {}'.format(error))

if __name__ == '__main__':
    main()


