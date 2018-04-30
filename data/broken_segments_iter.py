import numpy
import vtk

class BrokenSegmentsIter:

    def __init__(self, locator, brokenLine):
        """
        Constructor
        @param locator vtkCellLocator instance attached to the above grid
        @param brokenLine instance of BrokenLineIter
        """

        self.locator = locator
        self.data = []
        brokenLine.reset()
        for bl in brokenLine:
            t0, t1 = bl.getBegParamCoord(), bl.getEndParamCoord()
            dt = t1 - t0
            p0, p1 = bl.getBegPoint(), bl.getEndPoint()
            res = self.__collectLineGridIntersectionPoints(p0, p1, tol=1.e-10)
            for i in range(len(res) - 1):
                cellId, xia, ta = res[i]
                cellId, xib, tb = res[i + 1]
                self.data.append( (cellId, xia, xib, ta, tb) )

        self.numSegs = len(self.data)
        self.reset()


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
        if self.index < self.numSegs:
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

        eps = 1.234e-10

        # things we need to define
        cellId = vtk.mutable(-1)
        subId = vtk.mutable(-1)
        xi = numpy.zeros((3,), numpy.float64)
        tbar = vtk.mutable(-1.)
        cell = vtk.vtkGenericCell()
        pBeg3d = numpy.zeros((3,), numpy.float64)
        pEnd3d = numpy.zeros((3,), numpy.float64)
        point = numpy.zeros((3,), numpy.float64)
        weights = numpy.array((4,), numpy.float64)

        # VTK wants 3d positions
        # perturb the position to avoid muyltiple cells
        # to claim the same intersection point
        print p0
        pBeg3d[:2] = p0 - 0.67634534*eps
        pEnd3d[:2] = p1 + 0.48764787*eps

        res = []

        # add starting point
        cId = self.locator.FindCell(pBeg3d, tol, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2], 0.) )
        else:
            print('Warning: starting point {} not found!'.format(p0))
        tLast = 0.0
        pBeg3d += eps

        # find all intersection points in between
        found = True
        while found:

            found = self.locator.IntersectWithLine(pBeg3d, pEnd3d, tol, tbar, 
                                                   point, xi, subId, cellId)
            if found and tLast < 1.0 - eps:
                # correct the line param coord for the fact that we
                # moved the starting point
                t = tLast + (tbar.get() - tLast)/(1.0 - tLast)
                # store
                res.append( (cellId.get(), xi[:2], t) )
                # store the last line param coord
                tLast = t
                # reset the starting point of the ray
                pBeg3d[:2] = point[:2] + eps
            else:
                found = False
            
        # add last point 
        cId = self.locator.FindCell(pEnd3d, tol, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2], 1.) )
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
    parser.add_argument('-p', dest='points', default='', nargs='?', help='Points describing broken line in the format "(x0, y0),(x1, y1)..."')
    args = parser.parse_args()

    csr = CubedsphereReader(filename=args.input)
    points = eval(args.points)
    bl = BrokenLineIter(points)
    bs = BrokenSegmentsIter(csr.getUnstructuredGridCellLocator(), bl)
    for s in bs:
        print s.getCellId(), s.getBegCellParamCoord(), s.getEndCellParamCoord(), s.getBegLineParamCoord(), s.getEndLineParamCoord()

    

if __name__ == '__main__':
    main()
