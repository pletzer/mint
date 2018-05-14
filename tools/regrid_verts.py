from ugrid_reader import UgridReader
from latlon_reader import LatLonReader
from regrid_base import RegridBase
import numpy 
import vtk


class RegridVerts(RegridBase):

    # tolerance for finding cell
    EPS = 1.e-12

    def __init__(self):
        """
        Constructor
        no args
        """
        super(RegridVerts, self).__init__()


    def computeWeights(self):
        """
        Compute the interpolation weights
        assuming the field values are stored on vertices
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
                        self.weights[k] = self.ZERO4x4.copy()
                    # a point can only be in one cell, so no need to use +=
                    self.weights[k][i0, :] = ws


###############################################################################
def main():
    from math import pi, sin, cos, log, exp
    import argparse

    parser = argparse.ArgumentParser(description='Regriod edge field as if it were a nodal field')
    parser.add_argument('-s', dest='src', default='src.vtk', help='Specify source file in VTK unstructured grid format')
    parser.add_argument('-v', dest='varname', default='edge_integrated_velocity', help='Specify edge staggered field variable name in source VTK file')
    parser.add_argument('-d', dest='dst', default='dst.vtk', help='Specify destination file in VTK unstructured grid format')
    parser.add_argument('-o', dest='output', default='', help='Specify output VTK file where regridded edge data is saved')
    args = parser.parse_args()

    rgrd = RegridVerts()
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


