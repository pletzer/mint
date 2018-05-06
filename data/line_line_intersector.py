import numpy

class LineLineIntersector:

    def __init__(self):
        """
        Constructor
        no args
        """
        self.eps = 1.234e-12
        self.mat = numpy.zeros((2,2), numpy.float64)
        self.invMatTimesDet = numpy.zeros((2,2), numpy.float64)
        self.invMatTimesDetDotRhs = numpy.zeros((2,), numpy.float64)
        self.rhs = numpy.zeros((2,), numpy.float64)


    def setPoints(self, p0, p1, q0, q1):
        """
        Set points 
        @param p0 starting point of first line
        @param p1 end point of first line
        @param q0 starting point of second line
        @param q1 end point of second line
        """
        self.rhs[:] = q0 - p0
        self.mat[:, 0] = p1 - p0
        self.mat[:, 1] = q0 - q1
        self.invMatTimesDet[0, 0] = self.mat[1, 1]
        self.invMatTimesDet[1, 1] = self.mat[0, 0]
        self.invMatTimesDet[0, 1] = -self.mat[0, 1]
        self.invMatTimesDet[1, 0] = -self.mat[1, 0]
        self.solTimesDet = self.invMatTimesDet.dot(self.rhs)
        self.det = self.mat[0, 0]*self.mat[1, 1] - self.mat[0, 1]*self.mat[1, 0]

    def getDet(self):
        """
        Get the determinant
        @return determinant
        """
        return self.det

    def isSingular(self, tol):
        """
        Check if there is a solution
        @param tol tolerance
        @return True if there is one or more solutions
        """
        if abs(self.det) < tol:
            return True
        return False


    def hasSolution(self, tol):
        """
        Check if there is a solution
        @param tol tolerance
        @return True if there is one or more solutions
        """
        if abs(self.getDet()) > tol or self.solTimesDet.dot(self.solTimesDet) < tol:
            return True
        return False


    def solve(self):
        """
        Solve the system
        @param p0 starting point of first line
        @param p1 end point of first line
        @param q0 starting point of second line
        @param q1 end point of second line
        @return solution
        """
        return self.solTimesDet / self.det

###############################################################################
def test1():
    # standard
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([1., -1.])
    q1 = numpy.array([1., 2.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    xi1, xi2 = lli.solve()
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
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))

def test3():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([0., 1.])
    q1 = numpy.array([1., 1.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(not lli.hasSolution(1.e-10))

if __name__ == '__main__':
    test1()
    test2()
    test3()
