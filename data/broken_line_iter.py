import numpy

class BrokenLineIter:

    def __init__(self, points):
        """
        Constructor
        @param points list of points
        """
        self.index = 0
        self.numSegs = len(points) - 1

        lengths = numpy.array([numpy.linalg.norm(points[i + 1] - points[i]) \
                        for i in range(self.numSegs)])
        self.ts = numpy.array([0.] + list(lengths.cumsum()))
        self.ts /= self.ts[-1]
        print self.ts


    def __iter__(self):
        return self

    def next(self):
        """
        Update iterator
        """
        if self.index < self.numSegs:
            index = self.index
            self.index += 1
            return self
        else:
            raise StopIteration()

    def getBegParamCoord(self):
        """
        Get the parametric coordinate at the beginning of the segment
        @return t value
        """
        return self.ts[self.index - 1]
        

    def getEndParamCoord(self):
        """
        Get the parametric coordinate at the end of the segment
        @return t value
        """
        return self.ts[self.index]

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return self.index - 1


###############################################################################
def test():
    pts = numpy.array([[0., -1.], [0.5, 0.3], [-2., 1.1]])
    bl = BrokenLineIter(points=pts)
    for s in bl:
        print s.getIndex(), s.getBegParamCoord(), s.getEndParamCoord()


if __name__ == '__main__':
    test()
