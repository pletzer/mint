#include <mntLineLineIntersector.h>
#include <cassert>

void test1() {
    // standard
    double tol = 1.e-10;
    double p0[] = {0., 0.};
    double p1[] = {2., 0.};
    double q0[] = {1., -1.};
    double q1[] = {1., 2.};
    LineLineIntersector lli;
    lli.setPoints(p0, p1, q0, q1);
    Vector<double> xi = lli.getSolution();
    assert(abs(xi[0] - 1./2.) < tol);
    assert(abs(xi[1] - 1./3.) < tol);
}

/*

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
    print pa, pb
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

*/

int main(int argc, char** argv) {
    test1();
/*
    test2()
    test3()
    testNoOverlap()
    testNoOverlap2()
    testPartialOverlap()
    testPartialOverlap2()
    testPartialOverlap3()
    testQInsideP()
    testPInsideQ()
    */
    return 0;
}
