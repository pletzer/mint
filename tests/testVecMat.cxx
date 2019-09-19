#include <mntVec2.h>
#include <mntVec3.h>
#include <mntVec4.h>
#include <mntVec9.h>

void testVec2() {
    Vec2 v1;
    v1(0) = 10;
    v1(1) = 20;
    Vec2 v2;
    v2[0] = 1;
    v2[1] = 2;
    Vec2 w1 = v1 + v2;
    Vec2 w2 = v1 - v2;
    Vec2 w3 = 1.0 - v1;
    Vec2 w4 = exp(log(sin(v1 * v2 + 2.3) + 1.01));
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