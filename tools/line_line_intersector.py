import numpy

class LineLineIntersector:

    def __init__(self):
        """
        Constructor
        no args
        """
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
        self.p0, self.p1, self.q0, self.q1 = p0, p1, q0, q1
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


    def computeBegEndParamCoords(self):
        """
        Compute the begin/end parametric coordinates
        """
        dp = self.p1 - self.p0
        dp2 = dp.dot(dp)
        # lambda @ q0
        lm0 = (self.q0 - self.p0).dot(dp)/dp2
        # lambda @ q1
        lm1 = (self.q1 - self.p0).dot(dp)/dp2

        self.lamBeg = min(max(lm0, 0.0), 1.0)
        self.lamEnd = max(min(lm1, 1.0), 0.0)
        if self.lamEnd < self.lamBeg:
            # switch
            lbeg, lend = self.lamEnd, self.lamBeg
            self.lamBeg, self.lamEnd = lbeg, lend


    def getBegEndPoints(self):
        """
        Get the begin/end points of overlap
        @return tuple
        """
        dp = self.p1 - self.p0
        return self.p0 + self.lamBeg*dp, \
               self.p0 + self.lamEnd*dp


    def getBegEndParamCoords(self):
        """
        Get the begin/end parametric coordinates of overlap
        @return tuple
        """
        return self.lamBeg, self.lamEnd


    def hasSolution(self, tol):
        """
        Check if there is a solution
        @param tol tolerance
        @return True if there is one or more solutions
        """
        if abs(self.getDet()) > tol:
            return True
        if abs(self.solTimesDet.dot(self.solTimesDet)) < tol:
            # determinant is zero, p1 - p0 and q1 - q0 are on
            # the same ray
            self.computeBegEndParamCoords()
            if abs(self.lamEnd - self.lamBeg) > tol:
                return True

        return False


    def getSolution(self):
        """
        Get the solution
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
    xi1, xi2 = lli.getSolution()
    print('xi1 = {} xi2 = {}'.format(xi1, xi2))
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

def testNoOverlap():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([numpy.pi/2., 0.])
    q0 = numpy.array([-2., 0.])
    q1 = numpy.array([-1., 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(not lli.hasSolution(1.e-10))

def testNoOverlap2():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([numpy.pi/2., 0.])
    p1 = numpy.array([0., 0.])
    q0 = numpy.array([-2., 0.])
    q1 = numpy.array([-1., 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(not lli.hasSolution(1.e-10))


def testPartialOverlap():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([numpy.pi/2., 0.])
    q0 = numpy.array([-2., 0.])
    q1 = numpy.array([0.5, 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))
    pa, pb = lli.getBegEndPoints()
    u = (p1 - p0)/numpy.sqrt((p1 - p0).dot(p1 - p0))
    assert(abs((pa - p0).dot(u)) < 1.e-10)
    assert(abs((pb - q1).dot(u)) < 1.e-10)

def testPartialOverlap2():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([numpy.pi/2., 0.])
    q0 = numpy.array([0.1, 0.])
    q1 = numpy.array([numpy.pi, 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))
    pa, pb = lli.getBegEndPoints()
    u = (p1 - p0)/numpy.sqrt((p1 - p0).dot(p1 - p0))
    assert(abs((pa - q0).dot(u)) < 1.e-10)
    assert(abs((pb - p1).dot(u)) < 1.e-10)

def testPartialOverlap3():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([numpy.pi/2., 0.])
    p1 = numpy.array([0., 0.])
    q0 = numpy.array([0.1, 0.])
    q1 = numpy.array([numpy.pi, 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))
    pa, pb = lli.getBegEndPoints()
    u = (p1 - p0)/numpy.sqrt((p1 - p0).dot(p1 - p0))
    print('pa = {} pb = {}'.format(pa, pb))
    assert(abs((pa - p0).dot(u)) < 1.e-10)
    assert(abs((pb - q0).dot(u)) < 1.e-10)


def testQInsideP():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([1., 0.])
    q0 = numpy.array([0.1, 0.])
    q1 = numpy.array([0.8, 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))
    pa, pb = lli.getBegEndPoints()
    #print pa, pb
    u = (p1 - p0)/numpy.sqrt((p1 - p0).dot(p1 - p0))
    assert(abs((pa - q0).dot(u)) < 1.e-10)
    assert(abs((pb - q1).dot(u)) < 1.e-10)

def testPInsideQ():
    # no solution
    tol = 1.e-10
    p0 = numpy.array([0.1, 0.])
    p1 = numpy.array([0.9, 0.])
    q0 = numpy.array([0., 0.])
    q1 = numpy.array([1., 0.])
    lli = LineLineIntersector()
    lli.setPoints(p0, p1, q0, q1)
    det = lli.getDet()
    assert(abs(det) < 1.e-10)
    assert(lli.hasSolution(1.e-10))
    pa, pb = lli.getBegEndPoints()
    #print pa, pb
    u = (p1 - p0)/numpy.sqrt((p1 - p0).dot(p1 - p0))
    assert(abs((pa - p0).dot(u)) < 1.e-10)
    assert(abs((pb - p1).dot(u)) < 1.e-10)



if __name__ == '__main__':
    test1()
    test2()
    test3()
    testNoOverlap()
    testNoOverlap2()
    testPartialOverlap()
    testPartialOverlap2()
    testPartialOverlap3()
    testQInsideP()
    testPInsideQ()
