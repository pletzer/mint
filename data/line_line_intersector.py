import numpy

class LineLineIntersector:

	def __init__(self):
		
		self.mat = numpy.zeros((2,2), numpy.float64)
		self.rhs = numpy.zeros((2,), numpy.float64)
		self.sol = []


	def reset(self):
		self.rhs *= 0.0


	def setLine1(self, p0, p1):
		self.rhs -= p0
		self.mat[:, 0] = p1 - p0


	def setLine2(self, q0, q1):
		self.rhs += q0
		self.mat[:, 1] = q0 - q1


	def solve(self):
		self.sol = numpy.linalg.solve(self.mat, self.rhs)


	def getParamCoords(self):
		return self.sol

###############################################################################
def test1():
    tol = 1.e-10
    p0 = numpy.array([0., 0.])
    p1 = numpy.array([2., 0.])
    q0 = numpy.array([1., -1.])
    q1 = numpy.array([1., 2.])
    lli = LineLineIntersector()
    lli.reset()
    lli.setLine1(p0, p1)
    lli.setLine2(q0, p0)
    lli.solve()
    xi, eta = lli.getParamCoords()
    assert(abs(xi - 1.) < tol)
    assert(abs(eta - 1./3.) < tol)

if __name__ == '__main__':
    test1()
