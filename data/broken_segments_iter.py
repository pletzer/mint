import numpy
import vtk
from line_line_intersector import LineLineIntersector


class BrokenSegmentsIter:

    def __init__(self, grid, locator, brokenLine):
        """
        Constructor
        @param grid instance of vtkUnstructuredGrid
        @param locator vtkCellLocator instance attached to the above grid
        @param brokenLine instance of BrokenLineIter
        """

        # small tolerances 
        self.eps = 1.73654365e-12
        self.eps100 = 100. * self.eps
        self.tol = 1.e-2 # to determine if a point is inside a cell


        self.grid = grid
        self.locator = locator

        brokenLine.reset()
        for bl in brokenLine:
            t0, t1 = bl.getBegParamCoord(), bl.getEndParamCoord()
            dt = t1 - t0
            p0, p1 = bl.getBegPoint(), bl.getEndPoint()

            # res is  [ (cellId, xi, t), ...]
            res = self.__collectLineGridSegments(p0, p1)

            # DEBUG
            diffBeg = numpy.array([3.9269908169872414, 1.04089353704597, 0.]) - p0
            diffEnd = numpy.array([3.141592653589793, 1.1780972450961724, 0.]) - p1
            verbose = False
            if abs(diffBeg.dot(diffBeg)) < 1.e-8 and abs(diffEnd.dot(diffEnd)) < 1.e-8:
                verbose = True

            # re-arrange the data cellId -> [[t0, xi0], [t1, xi1], ...]
            c2s = {}
            for e in res:
                cId, xi, t = e
                c2s[cId] = c2s.get(cId, []) + [(t, xi)]

            # {(ta, tb) : (cellId, xia, xib, coeff), ...}
            data = {}
            for cId, v in c2s.items():
                v.sort(lambda x, y: cmp(x[0], y[0]))
                n = len(v)
                for i in range(n - 1):
                    ta, xia = v[i]
                    tb, xib = v[i + 1]
                    data[(ta, tb)] = [cId, xia, xib, 1.0]
            if verbose:
                print '%%%% data = ', data

        # turn data into a list [[(ta, tb), [cId, xia, xib, coeff]],...]
        self.data = [[k, v] for k, v in data.items()]

        # sort 
        self.data.sort( lambda x, y: cmp(x[0], y[0]) )

        # assign coefficients that account for duplicity, ie segments 
        # that are shared between two cells
        self.__assignCoefficientsToSegments()

        self.numSegs = len(self.data)

        self.totalT = 0
        for tatb, cxiaxibc in self.data:
            ta, tb = tatb
            coeff = cxiaxibc[3]
            self.totalT += (tb - ta) * coeff

        # reset the iterator
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
        return self.segment[1][0]


    def getBegCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segment[1][1]
        

    def getEndCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segment[1][2]
 

    def getBegLineParamCoord(self):
        """
        Get the current line parametric coordinates at the beginning of segment
        @return 2d array
        """
        return self.segment[0][0]
        

    def getEndLineParamCoord(self):
        """
        Get the current line parametric coordinates at the end of segment
        @return 2d array
        """
        return self.segment[0][1]

    def getCoefficient(self):
        """
        Get the coefficient accounting for duplicates
        @return coefficient
        """
        return self.segment[1][3]
 

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return self.index


    def __assignCoefficientsToSegments(self):

        n = len(self.data)

        # remove zero length segments
        for i in range(n - 1, -1, -1):
            s = self.data[i]
            ta, tb = s[0]
            if abs(tb - ta) < self.eps100:
                del self.data[i]

        # reduce contribution for overlapping segmmnts. If two 
        # segments overlap then the coefficient of first segment
        # is set to 1.0 - overlap/(tb - ta). Assumes overlap 
        # can only happen for pairs of segment
        n = len(self.data)
        for i0 in range(n - 1):
            i1 = i0 + 1
            s0 = self.data[i0]
            s1 = self.data[i1]
            ta0, tb0 = s0[0]
            ta1, tb1 = s1[0]
            overlap = max(0., min(tb0, tb1) - max(ta1, ta0))
            s0[1][3] = 1.0 - overlap/(tb0 - ta0)


    def __collectIntersectionPoints(self, pBeg, pEnd):
        """
        Collect all the intersection points
        @param pBeg starting point
        @param pEnd end point
        @return [(cellId, lambda, point), ...]
        @note lambda is the linear parametric coordinate along the line
        """

        res = []

        intersector = LineLineIntersector()
        cellIds = vtk.vtkIdList()
        cellIds1 = vtk.vtkIdList()
        cellIds2 = vtk.vtkIdList()
        ptIds = vtk.vtkIdList()

        dp = pEnd - pBeg

        # find all the cells intersected by the line
        self.locator.FindCellsAlongLine(pBeg, pEnd, self.tol, cellIds)
        """
        # direction of the ray
        u = dp / numpy.sqrt(dp.dot(dp))
        # up
        z = numpy.array([0., 0., 1.])
        # direction perpendicular to the ray
        v = numpy.cross(z, u)
        pBeg1, pEnd1 = pBeg + self.tol*v, pEnd + self.tol*v
        pBeg2, pEnd2 = pBeg - self.tol*v, pEnd - self.tol*v
        self.locator.FindCellsAlongLine(pBeg1, pEnd1, self.tol, cellIds1)
        self.locator.FindCellsAlongLine(pBeg2, pEnd2, self.tol, cellIds2)
        # add the cell from cellIds1 and cellIds2 to cellids
        for i in range(cellIds1.GetNumberOfIds()):
            cellIds.InsertUniqueId(cellIds1.GetId(i))
        for i in range(cellIds2.GetNumberOfIds()):
            cellIds.InsertUniqueId(cellIds2.GetId(i))
        """

        # collect the intersection points in between
        for i in range(cellIds.GetNumberOfIds()):

            cId = cellIds.GetId(i)

            self.grid.GetCellPoints(cId, ptIds)

            # iterate over the quads' edges
            for j0 in range(4):

                j1 = (j0 + 1) % 4
                v0 = numpy.array(self.grid.GetPoint(ptIds.GetId(j0)))
                v1 = numpy.array(self.grid.GetPoint(ptIds.GetId(j1)))

                # look for an intersection
                intersector.setPoints(pBeg[:2], pEnd[:2], v0[:2], v1[:2])
                if not intersector.hasSolution(self.eps):
                    continue

                if abs(intersector.getDet()) > self.eps:
                    # normal intersection, 1 solution
                    lambRay, lambEdg = intersector.getSolution()

                    # is it valid? Intersection must be within (p0, p1) and (q0, q1)
                    if lambRay >= 0. - self.eps100 and lambRay <= 1. + self.eps100 and \
                        lambEdg >= 0. - self.eps100 and lambEdg <= 1. + self.eps100:

                        point = pBeg + lambRay*dp
                        res.append( (cId, lambRay, point) )

                else:
                    # det is almost zero
                    # looks like the two lines (p0, p1) and (q0, q1) are overlapping
                    # add the starting/ending points 
                    lama, lamb = intersector.getBegEndParamCoords()
                    pa = pBeg + lama*dp
                    pb = pBeg + lamb*dp
                    res.append( (cId, lama, pa) )
                    res.append( (cId, lamb, pb) )

        return res


    def __collectLineGridSegments(self, p0, p1):
        """
        Collect all the line-grid intersection points
        @param p0 starting point of the line
        @param p1 end point of the line 
        @return list of [ (cellId, xi, t), ...]
        """

        res = []

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

        # VTK wants 3d positions
        pBeg[:] = p0[0], p0[1], 0.0
        pEnd[:] = p1[0], p1[1], 0.0

        # perturb the position to avoid a singular
        # system when looking for edge-line 
        # intersections
        """
        pBeg[0] += -self.eps * 1.86512432134
        pBeg[1] += +self.eps * 2.76354653243
        pEnd[0] += +self.eps * 1.96524543545
        pEnd[1] += -self.eps * 0.82875646565
        """

        deltaPos = pEnd - pBeg

        # add starting point
        cId = self.locator.FindCell(pBeg, self.eps, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 0.) )
        else:
            pass
            #print('Warning: starting point {} not found!'.format(p0))

        tLast = 0.0
        cIdLast = cId
        fullSegment = False

        # find all intersection points in between
        #pBeg += eps*deltaPos
        #pEnd -= eps*deltaPos

        # let's be generous with the collection of cells
        intersections = self.__collectIntersectionPoints(pBeg, pEnd)

        # find the cell id of the neighbouring cells
        for cId, lambRay, point in intersections:

            found = self.grid.GetCell(cId).EvaluatePosition(point, closestPoint, subId, xi, dist, weights)
            if found:
                res.append( (cId, xi[:2].copy(), lambRay) )
            else:
                print('Warning: param coord search failed point {} in cell {}'.format(point, cId))

            
        # add last point 
        cId = self.locator.FindCell(pEnd, self.eps, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 1.) )
        else:
            pass
            #print('Warning: end point {} not found!'.format(p1))

        return res

#################################################################################################

def main():
    import argparse
    from math import pi
    from ugrid_reader import UgridReader
    from broken_line_iter import BrokenLineIter

    parser = argparse.ArgumentParser(description='Break line into segments')
    parser.add_argument('-i', dest='input', default='mesh_C4.nc', help='Specify input file')
    parser.add_argument('-p', dest='points', default='(0., 0.),(2*pi,0.)', help='Points describing broken line as "(x0, y0),(x1, y1)..."')
    args = parser.parse_args()

    ur = UgridReader(filename=args.input)
    points = eval(args.points)
    bl = BrokenLineIter(points)
    
    count = 0
    for b in bl:
        t0 = b.getBegParamCoord()
        t1 = b.getEndParamCoord()
        print('line {} t = {} -> {}'.format(count, t0, t1))
        count += 1

    bs = BrokenSegmentsIter(ur.getUnstructuredGrid(), 
                            ur.getUnstructuredGridCellLocator(), bl)
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
