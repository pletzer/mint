import numpy

class LineLineIntersector:

    def __init__(self):
        """
        Constructor
        no args
        """
        
        self.mat = numpy.zeros((2,2), numpy.float64)
        self.rhs = numpy.zeros((2,), numpy.float64)
        self.sol = numpy.zeros((2,), numpy.float64)


    def reset(self):
        """
        Rest the matrix system
        """
        self.rhs *= 0.0


    def setLine1(self, p0, p1):
        """
        Set the first line
        @param p0 starting point
        @param p1 end point
        """
        self.rhs -= p0
        self.mat[:, 0] = p1 - p0


    def setLine2(self, q0, q1):
        """
        Set the second line
        @param q0 starting point
        @param q1 end point
        """
        self.rhs += q0
        self.mat[:, 1] = q0 - q1


    def solve(self):
        """
        Solve the system
        """
        self.sol = numpy.linalg.solve(self.mat, self.rhs)


    def getParamCoords(self):
        """
        Get the parametric coordinates of the solution vector
        @return vector
        """
        return self.sol

###############################################################################
def test1():
    # standard
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([1., -1.])
    q1 = numpy.array([1., 2.])
    lli = LineLineIntersector()
    lli.reset()
    lli.setLine1(p0, p1)
    lli.setLine2(q0, q1)
    lli.solve()
    xi1, xi2 = lli.getParamCoords()
    print xi1, xi2
    assert(abs(xi1 - 1./2.) < tol)
    assert(abs(xi2 - 1./3.) < tol)

def test2():
    # degenerate solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([0., 0.])
    q1 = numpy.array([1., 0.])
    lli = LineLineIntersector()
    lli.reset()
    lli.setLine1(p0, p1)
    lli.setLine2(q0, q1)
    try:
        lli.solve()
        xi1, xi2 = lli.getParamCoords()
    except:
        # error expected
        pass

def test3():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([0., 1.])
    q1 = numpy.array([1., 1.])
    lli = LineLineIntersector()
    lli.reset()
    lli.setLine1(p0, p1)
    lli.setLine2(q0, q1)
    try:
        lli.solve()
        xi1, xi2 = lli.getParamCoords()
    except:
        # error expected
        pass

if __name__ == '__main__':
    test1()
    test2()
    test3()
