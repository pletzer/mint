from regrid_base import RegridBase
import numpy 
import vtk


class RegridAvgVerts(RegridBase):

    EPS = 1.e-12

    def __init__(self):
        """
        Constructor
        no args
        """
        super(RegridAvgVerts, self).__init__()


    def computeWeights(self):
        """
        Compute the interpolation weights
        """

        dstPtIds = vtk.vtkIdList()
        cell = vtk.vtkGenericCell()
        pcoords = numpy.zeros((3,), numpy.float64)
        ws = numpy.zeros((4,), numpy.float64)
        
        numSrcCells = self.srcGrid.GetNumberOfCells()
        numDstCells = self.dstGrid.GetNumberOfCells()

        # iterate over the dst grid cells
        for dstCellId in range(numDstCells):

            # iterate over the four vertices of the dst cell
            self.dstGrid.GetCellPoints(dstCellId, dstPtIds)
            for i0 in range(4):

                dstVert = self.dstGrid.GetPoint(dstPtIds.GetId(i0))

                # bilinear interpolation
                srcCellId = self.srcLoc.FindCell(dstVert, self.EPS, cell, pcoords, ws)
                if srcCellId >= 0:
                    k = (dstCellId, srcCellId)
                    if not self.weights.has_key(k):
                        # initialize the weights
                        self.weights[k] = self.ZERO4x4
                    self.weights[k][i0, :] = ws


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

    parser = argparse.ArgumentParser(description='Regrid edge field by averaging vertex values on edges')
    parser.add_argument('-s', dest='src', default='mesh_C4.nc', help='Specify UGRID source grid file')
    parser.add_argument('-d', dest='dst', default='mesh_C4.nc', help='Specify UGRID destination grid file')
    parser.add_argument('-S', dest='srcStreamFunc', default='x', 
                        help='Stream function as a function of x (longitude in rad) and y (latitude in rad)')
    args = parser.parse_args()

    rgrd = RegridAvgVerts()
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


