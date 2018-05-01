import numpy
import vtk
from line_line_intersector import LineLineIntersector


class BrokenSegmentsIter:

    def __init__(self, grid, locator, brokenLine):
        """
        Constructor
        @param locator vtkCellLocator instance attached to the above grid
        @param brokenLine instance of BrokenLineIter
        """

        self.grid = grid
        self.locator = locator
        self.data = []
        brokenLine.reset()
        self.totalT = 0.0
        for bl in brokenLine:
            t0, t1 = bl.getBegParamCoord(), bl.getEndParamCoord()
            dt = t1 - t0
            p0, p1 = bl.getBegPoint(), bl.getEndPoint()
            res = self.__collectLineGridIntersectionPoints(p0, p1)
            # expect 2 or more points
            for i in range(len(res) - 1):
                cIda, xia, lama = res[i]
                cIdb, xib, lamb = res[i + 1]
                if cIda != cIdb:
                    print('Warning: cell {} != {} should not change between beg/end of segment'.format(cIda, cIdb))
                ta = t0 + lama*dt
                tb = t0 + lamb*dt
                self.totalT += tb - ta
                self.data.append( (cIda, xia, xib, ta, tb) )

        self.numSegs = len(self.data)
        self.reset()


    def getIntegratedParamCoord(self):
        """
        Get the integrated linear parametric coordinates
        @return value
        """
        return self.totalT


    def reset(self):
        """
        Reset the counter
        """
        self.index = -1
        self.segment = None


    def __iter__(self):
        return self


    def next(self):
        """
        Update iterator
        """
        if self.index < self.numSegs - 1:
            index = self.index
            self.index += 1
            self.segment = self.data[self.index]
            return self
        else:
            raise StopIteration()


    def getCellId(self):
        """
        Get the current cell Id
        @return index
        """
        return self.segment[0]


    def getBegCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segment[1]
        

    def getEndCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segment[2]
 

    def getBegLineParamCoord(self):
        """
        Get the current line parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segment[3]
        

    def getEndLineParamCoord(self):
        """
        Get the current line parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segment[4]
 

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return self.index


    def __collectLineGridIntersectionPoints(self, p0, p1, tol=1.e-10):
        """
        Collect all the line-grid intersection points
        @param p0 starting point of the line
        @param p1 end point of the line 
        @param tol tolerance
        @return list of [(cellId, xi, t), ...]
        """

        eps = 1.234e-12

        # things we need to define
        ptIds = vtk.vtkIdList()
        cell = vtk.vtkGenericCell()
        cellIds = vtk.vtkIdList()
        subId = vtk.mutable(-1)
        dist = vtk.mutable(0.)
        xi = numpy.zeros((3,), numpy.float64)
        pBeg = numpy.zeros((3,), numpy.float64)
        pEnd = numpy.zeros((3,), numpy.float64)
        point = numpy.zeros((3,), numpy.float64)
        closestPoint = numpy.zeros((3,), numpy.float64)
        weights = numpy.array((4,), numpy.float64)
        intersector = LineLineIntersector()

        # VTK wants 3d positions
        pBeg[:] = p0[0], p0[1], 0.0
        pEnd[:] = p1[0], p1[1], 0.0

        # perturb the position to avoid multiple cells
        # claiming the same intersection point
        pBeg[0] += eps*1.86512432134
        pBeg[1] -= eps*2.76354653243
        pEnd[0] -= eps*1.96524543545
        pEnd[1] += eps*0.82875646565

        deltaPos = pEnd - pBeg

        res = []

        # add starting point
        cId = self.locator.FindCell(pBeg, tol, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 0.) )
        else:
            print('Warning: starting point {} not found!'.format(p0))
        tLast = 0.0

        # find all intersection points in between
        pBeg += eps*deltaPos
        pEnd -= eps*deltaPos

        # find all the cells intersected by the line
        self.locator.FindCellsAlongLine(pBeg, pEnd, tol, cellIds)

        # iterate over the intersected cells
        for i in range(cellIds.GetNumberOfIds()):

            cId = cellIds.GetId(i)

            # get the point Ids of this cell
            self.grid.GetCellPoints(cId, ptIds)

            # iterate over the quads' edges
            for j in range(4):

                v0 = numpy.array(self.grid.GetPoint(ptIds.GetId(j)))
                v1 = numpy.array(self.grid.GetPoint(ptIds.GetId((j + 1) % 4)))

                # look for an intersection
                intersector.reset()
                intersector.setLine1(pBeg[:2], pEnd[:2])
                intersector.setLine2(v0[:2], v1[:2])
                intersector.solve()
                lambRay, lambEdg = intersector.getParamCoords()

                if lambRay >= 0. - eps and lambRay <= 1. + eps and \
                    lambEdg >= 0. - eps and lambEdg <= 1. + eps:

                    # found valid intersection, compute the cell 
                    # parametric coords
                    point = pBeg + lambRay*(pEnd - pBeg)
                    found = self.grid.GetCell(cId).EvaluatePosition(point, closestPoint, 
                                                                    subId, xi, 
                                                                    dist, weights)
                    if found:
                        # add
                        res.append( (cId, xi[:2].copy(), lambRay) )
                    else:
                        print('Warning: found intersection point {} but cell {} does not contain it'.format(point, cId))

            
        # add last point 
        cId = self.locator.FindCell(pEnd, tol, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 1.) )
        else:
            print('Warning: end point {} not found!'.format(p1))

        return res

def main():
    import argparse
    from math import pi
    from cubedsphere_reader import CubedsphereReader
    from broken_line_iter import BrokenLineIter

    parser = argparse.ArgumentParser(description='Break line into segments')
    parser.add_argument('-i', dest='input', default='mesh_C4.nc', help='Specify input file')
    parser.add_argument('-p', dest='points', default='', nargs='?', help='Points describing broken line as "(x0, y0),(x1, y1)..."')
    args = parser.parse_args()

    csr = CubedsphereReader(filename=args.input)
    points = eval(args.points)
    bl = BrokenLineIter(points)
    
    count = 0
    for b in bl:
        t0 = b.getBegParamCoord()
        t1 = b.getEndParamCoord()
        print('line {} t = {} -> {}'.format(count, t0, t1))
        count += 1

    bs = BrokenSegmentsIter(csr.getUnstructuredGrid(), csr.getUnstructuredGridCellLocator(), bl)
    count = 0
    for s in bs:
        cellId = s.getCellId()
        xia = s.getBegCellParamCoord()
        xib = s.getEndCellParamCoord()
        ta = s.getBegLineParamCoord()
        tb = s.getEndLineParamCoord()
        print('seg {} in cell {} t = {} -> {} xi = {} -> {}'.format(count, cellId, ta, tb, xia, xib))
        count += 1

    print('Integrated t = {}'.format(s.getIntegratedParamCoord()))

    

if __name__ == '__main__':
    main()
