from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
from regrid_base import RegridBase
from polyline_iter import PolylineIter
from polysegment_iter import PolysegmentIter
import numpy
import vtk


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
                id0 = dstPtIds.GetId(i0)
                id1 = dstPtIds.GetId(i1)
                dstEdgePt0 = self.dstGrid.GetPoint(id0)
                dstEdgePt1 = self.dstGrid.GetPoint(id1)

                # compute the intersections of the dst edge with the src grid
                psi = PolysegmentIter(self.srcGrid, self.srcLoc, dstEdgePt0, dstEdgePt1)

                # compute the contributions to this edge
                for seg in psi:

                    srcCellId = seg.getCellId()
                    xia = seg.getBegCellParamCoord()
                    xib = seg.getEndCellParamCoord()
                    coeff = seg.getCoefficient()

                    dxi = xib - xia
                    xiMid = 0.5*(xia + xib)

                    k = (dstCellId, srcCellId)

                    # compute the weights from each src edge
                    ws = numpy.array([+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                                      + dxi[1] * (0.0 + xiMid[0]) * coeff,
                                      + dxi[0] * (0.0 + xiMid[1]) * coeff,
                                      + dxi[1] * (1.0 - xiMid[0]) * coeff])

                    if not k in self.weights:
                        # initialize the weights
                        self.weights[k] = self.ZERO4x4.copy()

                    self.weights[k][i0, :] += ws

                totalT = psi.getIntegratedParamCoord()
                if abs(totalT - 1.0) > 1.e-6:
                    print('Warning: total t of segment: {:.3f} != 1 (diff={:.1g}), dst cell {} points=[{}, {}]'.format(\
                        totalT, totalT - 1.0, dstCellId, dstEdgePt0, dstEdgePt1))


###############################################################################
def main():
    from math import pi, sin, cos, log, exp
    import argparse

    parser = argparse.ArgumentParser(description='Regriod edge field')
    parser.add_argument('-s', dest='src', default='src.vtk', help='Specify source file in VTK unstructured grid format')
    parser.add_argument('-v', dest='varname', default='edge_integrated_velocity', help='Specify edge staggered field variable name in source VTK file')
    parser.add_argument('-d', dest='dst', default='dst.vtk', help='Specify destination file in VTK unstructured grid format')
    parser.add_argument('-o', dest='output', default='', help='Specify output VTK file where regridded edge data is saved')
    args = parser.parse_args()

    rgrd = RegridEdges()
    rgrd.setSrcFile(args.src)
    rgrd.setDstFile(args.dst)
    rgrd.computeWeights()

    # compute edge integrals
    srcEdgeVel = rgrd.getSrcEdgeData(args.varname)

    # regrid/apply the weights 
    dstEdgeVel = rgrd.applyWeights(srcEdgeVel)

    # loop integrals for each cell
    cellLoops = dstEdgeVel.sum(axis=1)

    # statistics
    absCellIntegrals = numpy.abs(cellLoops)
    minAbsLoop = absCellIntegrals.min()
    maxAbsLoop = absCellIntegrals.max()
    avgAbsLoop = absCellIntegrals.sum() / float(absCellIntegrals.shape[0])

    print('Min/avg/max cell loop integrals: {}/{}/{}'.format(minAbsLoop, avgAbsLoop, maxAbsLoop))
    if args.output:
        rgrd.saveDstEdgeData(args.varname, dstEdgeVel, args.output)

if __name__ == '__main__':
    main()


