#include <mntLineTriangleIntersector.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

void test1() {
    // standard
    double tol = 1.e-10;

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
    std::cout << "test1: det = " << det << '\n';

    Vector<double> sol = lli.getSolution();
    double lam = sol[0];
    assert(abs(lam - 1.0) < tol);
    std::vector<double> xi = std::vector<double>(&sol[1], &sol[3]);
    std::cout << "test1: xi = " << xi[0] << ", " << xi[1] << '\n';
    assert(abs(xi[0] - 0.) < tol);
    assert(abs(xi[1] - 0.) < tol);
}


int main(int argc, char** argv) {

    test1();

    return 0;
}
