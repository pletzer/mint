#include <MvVector3.h>

int main() {
    
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

    return 0;
}