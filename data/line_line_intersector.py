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


    def solve(self, p0, p1, q0, q1):
        """
        Solve the system
        @param p0 starting point of first line
        @param p1 end point of first line
        @param q0 starting point of second line
        @param q1 end point of second line
        @return solution
        """
        self.rhs = q0 - p0
        self.mat[:, 0] = p1 - p0
        self.mat[:, 1] = q0 - q1
        return numpy.linalg.solve(self.mat, self.rhs)

###############################################################################
def test1():
    # standard
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([1., -1.])
    q1 = numpy.array([1., 2.])
    lli = LineLineIntersector()
    xi1, xi2 = lli.solve(p0, p1, q0, q1)
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
    try:
        xi1, xi2 = lli.solve(p0, p1, q0, q1)
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
    try:
        xi1, xi2 = lli.solve(p0, p1, q0, q1)
    except:
        # error expected
        pass

if __name__ == '__main__':
    test1()
    test2()
    test3()
