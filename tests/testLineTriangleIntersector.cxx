#include <mntLineTriangleIntersector.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

#define TOL 1.e-10

void testNormal() {

    // line
    double p0[] = {0., 0., 1.};
    double p1[] = {0., 0., 0.};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lli;
    lli.setPoints(p0, p1, q0, q1, q2);

    double det = lli.getDet();
    assert(abs(det - (-1.)) < TOL);

    Vector<double> sol = lli.getSolution();
    double lam = sol[0];
    assert(abs(lam - 1.0) < TOL);
    std::vector<double> xi = std::vector<double>(&sol[1], &sol[3]);
    assert(abs(xi[0] - 0.) < TOL);
    assert(abs(xi[1] - 0.) < TOL);
}

void testNoSolution() {

    // line
    double p0[] = {0., 0., 0.1};
    double p1[] = {0., 1., 0.1};

    // triangle
    double q0[] = {0., 0., 0.};
    double q1[] = {1., 0., 0.};
    double q2[] = {0., 1., 0.};
    LineTriangleIntersector lli;
    lli.setPoints(p0, p1, q0, q1, q2);

    double det = lli.getDet();
    assert(abs(det - (0.)) < TOL);
    assert(!lli.hasSolution(TOL));
}



int main(int argc, char** argv) {

    testNormal();
    testNoSolution();

    return 0;
}
