#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>

#include <mntVec2.h>
#include <mntVec3.h>
#include <mntVec4.h>
#include <mntVec9.h>
#include <mntMat2x2.h>
#include <mntMat3x3.h>

const double EPS = 1.e-10;

void testVec2() {
    Vec2 v1;
    v1(0) = 10;
    v1(1) = 20;
    Vec2 v2;
    v2[0] = 1;
    v2[1] = 2;

    // test vector operations
    Vec2 w1 = v1 + v2;
    Vec2 w2 = v1 - v2;
    Vec2 w3 = 1.0 - v1;
    Vec2 w4 = exp(log(sin(v1 * v2 + 2.3) + 1.01));
    for (size_t i = 0; i < w1.size(); ++i) {
        assert(fabs(w1[i] - (v1[i] + v2[i])) < EPS);
        assert(fabs(w2[i] - (v1[i] - v2[i])) < EPS);
        assert(fabs(w3[i] - (1.0 - v1[i])) < EPS);
    }

    // test dot product of vectors
    double d = dot(v1, v2);
    assert(fabs(d - (v1[0]*v2[0] + v1[1]*v2[1])) < EPS);

    // test matrix dot operations
    Mat2x2 m(0);
    m(0, 0) = -1; m(0, 1) = 1;
    m(1, 0) = 2; m(1, 1) = 3;
    Vec2 x = dot(v1, m);
    Vec2 y = dot(m, v1);
    Mat2x2 m3 = dot(m, m);
    assert(fabs(x[0] - (v1[0]*m(0, 0) + v1[1]*m(1, 0))) < EPS);
    assert(fabs(x[1] - (v1[0]*m(0, 1) + v1[1]*m(1, 1))) < EPS);
    assert(fabs(y[0] - (m(0, 0)*v1[0] + m(0, 1)*v1[1])) < EPS);
    assert(fabs(y[1] - (m(1, 0)*v1[0] + m(1, 1)*v1[1])) < EPS);

    // test transpose
    Mat2x2 m4 = transpose(m);
    assert(fabs(m4(0, 0) - m(0, 0)) < EPS);
    assert(fabs(m4(0, 1) - m(1, 0)) < EPS);
    assert(fabs(m4(1, 0) - m(0, 1)) < EPS);
    assert(fabs(m4(1, 1) - m(1, 1)) < EPS);
}

void testVec3() {
    Vec3 v1;
    v1(0) = 10;
    v1(1) = 20;
    v1(2) = 30;
    Vec3 v2;
    v2[0] = 1;
    v2[1] = 2;
    v2[2] = 3;
    Vec3 w1 = v1 + v2;
    Vec3 w2 = v1 - v2;
    Vec3 w3 = 1.0 - v1;
    Vec3 w4 = exp(log(sin(v1 * v2 + 2.3) + 1.01));

    // test matrix dot operations
    Mat3x3 m(0);
    m(0, 0) = -1; m(0, 1) = 2; m(0, 2) = 3.;
    m(1, 0) = -4; m(1, 1) = 5; m(1, 2) = 6.;
    m(2, 0) = -7; m(2, 1) = 8; m(2, 2) = 9.;
    Vec3 x = dot(v1, m);
    Vec3 y = dot(m, v1);
    Mat3x3 m3 = dot(m, m);
    assert(fabs(x[0] - (v1[0]*m(0, 0) + v1[1]*m(1, 0) + v1[2]*m(2, 0))) < EPS);
    assert(fabs(x[1] - (v1[0]*m(0, 1) + v1[1]*m(1, 1) + v1[2]*m(2, 1))) < EPS);
    assert(fabs(x[2] - (v1[0]*m(0, 2) + v1[1]*m(1, 2) + v1[2]*m(2, 2))) < EPS);
    assert(fabs(y[0] - (m(0, 0)*v1[0] + m(0, 1)*v1[1] + m(0, 2)*v1[2])) < EPS);
    assert(fabs(y[1] - (m(1, 0)*v1[0] + m(1, 1)*v1[1] + m(1, 2)*v1[2])) < EPS);
    assert(fabs(y[2] - (m(2, 0)*v1[0] + m(2, 1)*v1[1] + m(2, 2)*v1[2])) < EPS);

    // test transpose
    Mat3x3 m4 = transpose(m);
    assert(fabs(m4(0, 0) - m(0, 0)) < EPS);
    assert(fabs(m4(0, 1) - m(1, 0)) < EPS);
    assert(fabs(m4(0, 2) - m(2, 0)) < EPS);
    assert(fabs(m4(1, 0) - m(0, 1)) < EPS);
    assert(fabs(m4(1, 1) - m(1, 1)) < EPS);
    assert(fabs(m4(1, 2) - m(2, 1)) < EPS);
    assert(fabs(m4(2, 0) - m(0, 2)) < EPS);
    assert(fabs(m4(2, 1) - m(1, 2)) < EPS);
    assert(fabs(m4(2, 2) - m(2, 2)) < EPS);
}

void testVec4() {
    Vec4 v1;
    v1(0) = 10;
    v1(1) = 20;
    v1(2) = 30;
    v1(3) = 40;
    Vec4 v2;
    v2[0] = 1;
    v2[1] = 2;
    v2[2] = 3;
    v2[3] = 4;
    Vec4 w1 = v1 + v2;
    Vec4 w2 = v1 - v2;
    Vec4 w3 = 1.0 - v1;
    Vec4 w4 = exp(log(sin(v1 * v2 + 2.3) + 1.01));
}

void testVec9() {
    Vec9 v1;
    v1(0) = 10;
    v1(1) = 20;
    v1(2) = 30;
    v1(3) = 40;
    v1(4) = 50;
    v1(5) = 60;
    v1(6) = 70;
    v1(7) = 80;
    v1(8) = 90;
    Vec9 v2;
    v2[0] = 1;
    v2[1] = 2;
    v2[2] = 3;
    v2[3] = 4;
    v2[4] = 5;
    v2[5] = 6;
    v2[6] = 7;
    v2[7] = 8;
    v2[8] = 9;
    Vec9 w1 = v1 + v2;
    Vec9 w2 = v1 - v2;
    Vec9 w3 = 1.0 - v1;
    Vec9 w4 = exp(log(sin(v1 * v2 + 2.3) + 1.01));
}

int main() {
    
    testVec2();
    testVec3();
    testVec4();
    testVec9();

    return 0;
}