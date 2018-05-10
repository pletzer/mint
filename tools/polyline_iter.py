import numpy

class PolylineIter:

    def __init__(self, points):
        """
        Constructor
        @param points list of points
        """
        self.points = [numpy.array(p) for p in points]
        self.numSegs = len(points) - 1

        lengths = numpy.array([numpy.linalg.norm(self.points[i + 1] - self.points[i]) \
                        for i in range(self.numSegs)])
        self.ts = numpy.array([0.] + list(lengths.cumsum()))
        self.ts /= self.ts[-1]
        self.reset()


    def reset(self):
        """
        Reset the counter
        """
        self.index = -1


    def __iter__(self):
        return self

    def next(self):
        """
        Update iterator
        """
        if self.index < self.numSegs - 1:
            self.index += 1
            return self
        else:
            raise StopIteration()

    def getBegPoint(self):
        """
        Get the point at the beginning of the line
        @return point
        """
        return self.points[self.index]


    def getEndPoint(self):
        """
        Get the point at the end of the line
        @return point
        """
        return self.points[self.index + 1]


    def getBegParamCoord(self):
        """
        Get the parametric coordinate at the beginning of the segment
        @return t value
        """
        return self.ts[self.index]
        

    def getEndParamCoord(self):
        """
        Get the parametric coordinate at the end of the segment
        @return t value
        """
        return self.ts[self.index + 1]

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return self.index


###############################################################################
def test():
    pts = numpy.array([[0., -1.], [0.5, 0.3], [-2., 1.1]])
    bl = PolylineIter(points=pts)
    bl.reset()
    for s in bl:
        print s.getIndex(), s.getBegParamCoord(), s.getEndParamCoord()


if __name__ == '__main__':
    test()
