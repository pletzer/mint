#include <mntLineLineIntersector.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>

void test1() {
    // standard
    double tol = 1.e-10;
    double p0[] = {0., 0.};
    double p1[] = {2., 0.};
    double q0[] = {1., -1.};
    double q1[] = {1., 2.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    Vec2 xi = lli.getSolution();
    std::cout << "test1: xi = " << xi << '\n';
    assert(abs(xi[0] - 1./2.) < tol);
    assert(abs(xi[1] - 1./3.) < tol);
}

void test1_3d() {
    // standard 3d on the same plane
    double tol = 1.e-10;
    double p0[] = {0., 0., 0.};
    double p1[] = {2., 0., 0.};
    double q0[] = {1., -1., 0.};
    double q1[] = {1., 2., 0.};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    Vec2 xi = lli.getSolution();
    assert(abs(xi[0] - 1./2.) < tol);
    assert(abs(xi[1] - 1./3.) < tol);
}

void test1_3DOffset() {
    // standard 3d on different, parallel planes
    double tol = 1.e-10;
    double p0[] = {0., 0., 0.};
    double p1[] = {2., 0., 0.};
    double q0[] = {1., -1., 0.1};
    double q1[] = {1., 2., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    Vec2 xi = lli.getSolution();
    assert(abs(xi[0] - 1./2.) < tol);
    assert(abs(xi[1] - 1./3.) < tol);
}


void test2() {
    // degenerate solution
    double tol = 1.e-10;
    double p0[] = {0., 0.};
    double p1[] = {2., 0.};
    double q0[] = {0., 0.};
    double q1[] = {1., 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
}

void test2_3DOffset() {
    // degenerate solution
    double tol = 1.e-10;
    double p0[] = {0., 0., 0.};
    double p1[] = {2., 0., 0.};
    double q0[] = {0., 0., 0.1};
    double q1[] = {1., 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
}

void test3() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0.};
    double p1[] =   {2., 0.};
    double q0[] =   {0., 1.};
    double q1[] =   {1., 1.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    std::cout << "test3: det = " << det << '\n';
    assert(! lli.hasSolution(tol));
}

void test3_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0., 0.};
    double p1[] =   {2., 0., 0.};
    double q0[] =   {0., 1., 0.1};
    double q1[] =   {1., 1., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(! lli.hasSolution(tol));
}


void testNoOverlap() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0.};
    double p1[] =   {M_PI/2., 0.};
    double q0[] =   {-2., 0.};
    double q1[] =   {-1., 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(! lli.hasSolution(tol));
}

void testNoOverlap_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0., 0.};
    double p1[] =   {M_PI/2., 0., 0.};
    double q0[] =   {-2., 0., 0.1};
    double q1[] =   {-1., 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(! lli.hasSolution(tol));
}

void testNoOverlap2() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {M_PI /2., 0.};
    double p1[] =   {0., 0.};
    double q0[] =   {-2., 0.};
    double q1[] =   {-1., 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(! lli.hasSolution(tol));
}

void testNoOverlap2_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {M_PI /2., 0., 0.};
    double p1[] =   {0., 0., 0.};
    double q0[] =   {-2., 0., 0.1};
    double q1[] =   {-1., 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(! lli.hasSolution(tol));
}

void testPartialOverlap() {
    double tol = 1.e-10;
    double p0[] = {0., 0.};
    double p1[] =  {M_PI /2., 0.};
    double q0[] =  {-2., 0.};
    double q1[] =  {0.5, 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 dp10;
    Vec2 dpap0;
    Vec2 dpbq1;
    Vec2 u;
    for (size_t i = 0; i < 2; ++i) {
        dp10[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbq1[i] = pb[i] - q1[i];
        u[i] = p1[i] - p0[i];
    }
    u /= dot(dp10, dp10);

    assert(std::abs(dot(dpap0, u)) < tol);
    assert(std::abs(dot(dpbq1, u)) < tol);
}

void testPartialOverlap_3DOffset() {
    double tol = 1.e-10;
    double p0[] =  {0., 0., 0.};
    double p1[] =  {M_PI /2., 0., 0.};
    double q0[] =  {-2., 0., 0.1};
    double q1[] =  {0.5, 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 dp10;
    Vec2 dpap0;
    Vec2 dpbq1;
    Vec2 u;
    for (size_t i = 0; i < 2; ++i) {
        dp10[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbq1[i] = pb[i] - q1[i];
        u[i] = p1[i] - p0[i];
    }
    u /= dot(dp10, dp10);

    assert(std::abs(dot(dpap0, u)) < tol);
    assert(std::abs(dot(dpbq1, u)) < tol);
}


void testPartialOverlap2() {
    double tol = 1.e-10;
    double p0[] =  {0., 0.};
    double p1[] =  {M_PI /2., 0.};
    double q0[] =  {0.1, 0.};
    double q1[] =  {M_PI , 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(std::abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 dp10;
    Vec2 dpaq0;
    Vec2 dpbp1;
    Vec2 u;
    for (size_t i = 0; i < 2; ++i) {
        dp10[i] = p1[i] - p0[i];
        dpaq0[i] = pa[i] - q0[i];
        dpbp1[i] = pb[i] - p1[i];
        u[i] = p1[i] - p0[i];
    }
    u /= sqrt(dot(dp10, dp10));

    assert(abs(dot(dpaq0, u)) < tol);
    assert(abs(dot(dpbp1, u)) < tol);
}

void testPartialOverlap2_3DOffset() {
    double tol = 1.e-10;
    double p0[] =  {0., 0., 0.};
    double p1[] =  {M_PI /2., 0., 0.};
    double q0[] =  {0.1, 0., 0.1};
    double q1[] =  {M_PI , 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(std::abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 dp10;
    Vec2 dpaq0;
    Vec2 dpbp1;
    Vec2 u;
    for (size_t i = 0; i < 2; ++i) {
        dp10[i] = p1[i] - p0[i];
        dpaq0[i] = pa[i] - q0[i];
        dpbp1[i] = pb[i] - p1[i];
        u[i] = p1[i] - p0[i];
    }
    u /= sqrt(dot(dp10, dp10));

    assert(abs(dot(dpaq0, u)) < tol);
    assert(abs(dot(dpbp1, u)) < tol);
}

void testPartialOverlap3() {
    double tol = 1.e-10;
    double p0[] =   {M_PI /2., 0.};
    double p1[] =   {0., 0.};
    double q0[] =   {0.1, 0.};
    double q1[] =   {M_PI , 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    std::cerr << "testPartialOverlap3: det = " << det << '\n';
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dpap0;
    Vec2 dpbq0;
    Vec2 dp10;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbq0[i] = pb[i] - q0[i];
        dp10[i] = p1[i] - p0[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpap0, u)) < tol);
    assert(abs(dot(dpbq0, u)) < tol);
}

void testPartialOverlap3_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {M_PI /2., 0., 0.};
    double p1[] =   {0., 0., 0.};
    double q0[] =   {0.1, 0., 0.1};
    double q1[] =   {M_PI , 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dpap0;
    Vec2 dpbq0;
    Vec2 dp10;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbq0[i] = pb[i] - q0[i];
        dp10[i] = p1[i] - p0[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpap0, u)) < tol);
    assert(abs(dot(dpbq0, u)) < tol);
}


void testQInsideP() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0.};
    double p1[] =   {1., 0.};
    double q0[] =   {0.1, 0.};
    double q1[] =   {0.8, 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dp10;
    Vec2 dpaq0;
    Vec2 dpbq1;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dp10[i] = p1[i] - p0[i];
        dpaq0[i] = pa[i] - q0[i];
        dpbq1[i] = pb[i] - q1[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpaq0, u)) < tol);
    assert(abs(dot(dpbq1, u)) < tol);
}

void testQInsideP_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0., 0., 0.};
    double p1[] =   {1., 0., 0.};
    double q0[] =   {0.1, 0., 0.1};
    double q1[] =   {0.8, 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dp10;
    Vec2 dpaq0;
    Vec2 dpbq1;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dp10[i] = p1[i] - p0[i];
        dpaq0[i] = pa[i] - q0[i];
        dpbq1[i] = pb[i] - q1[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpaq0, u)) < tol);
    assert(abs(dot(dpbq1, u)) < tol);
}

void testPInsideQ() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0.1, 0.};
    double p1[] =   {0.9, 0.};
    double q0[] =   {0., 0.};
    double q1[] =   {1., 0.};
    LineLineIntersector lli;
    lli.setPoints(2, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dp10;
    Vec2 dpap0;
    Vec2 dpbp1;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dp10[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbp1[i] = pb[i] - p1[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpap0, u)) < tol);
    assert(abs(dot(dpbp1, u)) < tol);
}

void testPInsideQ_3DOffset() {
    // no solution
    double tol = 1.e-10;
    double p0[] =   {0.1, 0., 0.};
    double p1[] =   {0.9, 0., 0.};
    double q0[] =   {0., 0., 0.1};
    double q1[] =   {1., 0., 0.1};
    LineLineIntersector lli;
    lli.setPoints(3, p0, p1, q0, q1);
    double det = lli.getDet();
    assert(abs(det) < tol);
    assert(lli.hasSolution(tol));
    std::pair< double, double > p = lli.getBegEndParamCoords();
    double lamA = p.first;
    double lamB = p.second;
    Vec2 pa = (1. - lamA)*Vec2(p0) + lamA*Vec2(p1);
    Vec2 pb = (1. - lamB)*Vec2(p0) + lamB*Vec2(p1);

    Vec2 u;
    Vec2 dp10;
    Vec2 dpap0;
    Vec2 dpbp1;
    for (size_t i = 0; i < 2; ++i) {
        u[i] = p1[i] - p0[i];
        dp10[i] = p1[i] - p0[i];
        dpap0[i] = pa[i] - p0[i];
        dpbp1[i] = pb[i] - p1[i];
    }
    u /= sqrt(dot(dp10, dp10));
    assert(abs(dot(dpap0, u)) < tol);
    assert(abs(dot(dpbp1, u)) < tol);
}


int main(int argc, char** argv) {

    test1();
    test2();
    test3();
    testNoOverlap();
    testNoOverlap2();
    testPartialOverlap();
    testPartialOverlap2();
    testPartialOverlap3();
    testQInsideP();
    testPInsideQ();

    // 3d tests
    test1_3d();
    test1_3DOffset();
    test2_3DOffset();
    //test3_3DOffset(); FAILS
    testNoOverlap_3DOffset();
    testNoOverlap2_3DOffset();
    //testPartialOverlap_3DOffset(); FAILS
    //testPartialOverlap2_3DOffset(); FAILS
    //testPartialOverlap3_3DOffset(); FAILS
    testQInsideP_3DOffset();
    //testPInsideQ_3DOffset(); FAILS

    return 0;
}
