#include <mntLineTriangleIntersector.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

#define TOL 1.e-10

void testGetParamLocation1() {

    Vec3 q0(3); q0[0] = 0.; q0[1] = 0.; q0[2] = 0.;
    Vec3 q1(3); q1[0] = 1.; q1[1] = 0.; q1[2] = 0.;
    Vec3 q2(3); q2[0] = 0.; q2[1] = 1.; q2[2] = 0.;

    Vec3 p(3); p[0] = 0.; p[1] = 0.1; p[2] = 0.;

    Vec2 xi = getTriangleParamLocation(q0, q1, q2, p);
    std::cerr << "testGetParamLocation1: xi = " << xi[0] << ',' << xi[1] << '\n';
    assert(std::abs(xi[0] - 0.) < TOL);
    assert(std::abs(xi[1] - 0.1) < TOL);
}

void testGetParamLocation2() {

    Vec3 q0(3); q0[0] = 0.; q0[1] = -90.; q0[2] = -90.;
    Vec3 q1(3); q1[0] = 0.; q1[1] = -90.; q1[2] = 0.;
    Vec3 q2(3); q2[0] = 0.; q2[1] = -45.; q2[2] = -90.;

    Vec3 p(3); p[0] = 0.; p[1] = -45.; p[2] = -90.;

    Vec2 xi = getTriangleParamLocation(q0, q1, q2, p);
    std::cerr << "testGetParamLocation2: xi = " << xi[0] << ',' << xi[1] << '\n';
    assert(std::abs(xi[0] - 0.) < TOL);
    assert(std::abs(xi[1] - 1.) < TOL);
}

void testGetParamLocation3() {

    Vec3 q0(3); q0[0] = 0.; q0[1] = -90.; q0[2] = -90.;
    Vec3 q1(3); q1[0] = 0.; q1[1] = -90.; q1[2] = 0.;
    Vec3 q2(3); q2[0] = 0.; q2[1] = -45.; q2[2] = -90.;

    Vec3 p(3); p[0] = 0.; p[1] = -90.; p[2] = -180.;

    Vec2 xi = getTriangleParamLocation(q0, q1, q2, p);
    std::cerr << "testGetParamLocation3: xi = " << xi[0] << ',' << xi[1] << '\n';
    assert(std::abs(xi[0] - (-1.)) < TOL);
    assert(std::abs(xi[1] - 0.) < TOL);
}

void testGeneral(const double pa[], const double pb[], 
                 const double q0[], const double q1[], const double q2[], 
                 double expectedLambda, const double expectedXi[]) {
    LineTriangleIntersector lti;
    lti.setPoints(pa, pb, q0, q1, q2);
    std::cerr << "testGeneral: det = " << lti.getDet() << '\n';
    assert(lti.hasSolution(TOL));
    Vec3 sol = lti.getSolution();
    assert(std::abs(sol[0] - expectedLambda) < TOL);
    assert(std::abs(sol[1] - expectedXi[0]) < TOL);
    assert(std::abs(sol[2] - expectedXi[1]) < TOL);
}

void testNormal() {

    // line
    double p0[] = {0., 0., 1.};
    double p1[] = {0., 0., 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (-1.)) < TOL);
    assert(lti.hasSolution(TOL));
    Vec3 sol = lti.getSolution();
    double lam = sol[0];
    assert(std::abs(lam - 1.0) < TOL);
    Vec2 xi = Vec2(&sol[1]);
    assert(std::abs(xi[0] - 0.) < TOL);
    assert(std::abs(xi[1] - 0.) < TOL);
}

void testNoSolution() {

    // line
    double p0[] = {0., 0., 0.1};
    double p1[] = {0., 1., 0.1};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    assert(!lti.hasSolution(TOL));
}

void testDegenerate1() {

    // line
    double p0[] = {0., 0.1, 0.};
    double p1[] = {1., 0.1, 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    assert(lti.hasSolution(TOL));
    const std::pair< double, double > lamBegEnd = lti.getBegEndParamCoords();
    assert(std::abs(lamBegEnd.first - 0.0) < TOL);
    assert(std::abs(lamBegEnd.second - 0.9) < TOL);    
}

void testDegenerate2() {

    // line
    double p0[] = {0.2, 0.1, 0.};
    double p1[] = {1.2, 0.1, 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    assert(lti.hasSolution(TOL));
    const std::pair< double, double > lamBegEnd = lti.getBegEndParamCoords();
    assert(std::abs(lamBegEnd.first - 0.0) < TOL);
    assert(std::abs(lamBegEnd.second - 0.7) < TOL);
}

void testDegenerate3() {

    // line
    double p0[] = {0., 0.1, 0.};
    double p1[] = {1.2, 0.1, 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    assert(lti.hasSolution(TOL));
    const std::pair< double, double > lamBegEnd = lti.getBegEndParamCoords();
    assert(std::abs(lamBegEnd.first - 0.0) < TOL);
    assert(std::abs(lamBegEnd.second - 0.75) < TOL);
}

void testDegenerate4() {

    // line
    double p0[] = {0., 0.1, 0.};
    double p1[] = {0.9, 0.1, 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    assert(lti.hasSolution(TOL));
    const std::pair< double, double > lamBegEnd = lti.getBegEndParamCoords();
    assert(std::abs(lamBegEnd.first - 0.0) < TOL);
    assert(std::abs(lamBegEnd.second - 1.0) < TOL);
}

void testCoPlanarNoSolution() {

    // line
    double p0[] = {1.2, 0.1, 0.};
    double p1[] = {2.0, 0.1, 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lti;
    lti.setPoints(p0, p1, q0, q1, q2);

    double det = lti.getDet();
    assert(std::abs(det - (0.)) < TOL);
    bool hasSol = lti.hasSolution(TOL);
    assert(!lti.hasSolution(TOL));
}

int main(int argc, char** argv) {


    testGetParamLocation3();
    testGetParamLocation2();
    testGetParamLocation1();

    double pa[] = {0., 0., 0.};
    double pb[] = {0., 0., 0.};
    double q0[] = {0., 0., 0.};
    double q1[] = {0., 0., 0.};
    double q2[] = {0., 0., 0.};
    double xi[] = {0., 0.};

    pa[0] = 0.; pa[1] = -90.; pa[2] = -180.;
    pb[0] = 0.; pb[1] = +90.; pb[2] = +180.;
    q0[0] = 0.; q0[1] = 0.; q0[2] = -60.;
    q1[0] = 0.; q1[1] = 0.; q1[2] = +60.;
    q2[0] = 1.; q2[1] = 0.; q2[2] = +60.;
    xi[0] = 0.5; xi[1] = 0.;
    testGeneral(pa, pb, q0, q1, q2, 0.5, xi);


    testNormal();
    testNoSolution();
    testDegenerate1();
    testDegenerate2();
    testDegenerate3();
    testDegenerate4();
    testCoPlanarNoSolution();

    return 0;
}
