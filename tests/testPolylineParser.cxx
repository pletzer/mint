#include "mntPolylineParser.h"
#undef NDEBUG // turn on asserts
#include <cassert>


void test2d() {
    PolylineParser pp(2);
    pp.parse("(1,2), (-1.2, 3.45), (-1.2e-3, 3.45e19)");
    const std::vector< Vec3 >& points = pp.getPoints();
    pp.print();
    assert(points.size() == 3);
}

void test3d() {
    PolylineParser pp(3);
    pp.parse("(1,2,3), (-1.2, 3.45, +6.79), (-1.2e-3, 3.45e19, -689.e-12)");
    const std::vector< Vec3 >& points = pp.getPoints();
    pp.print();
    assert(points.size() == 3);
}

void testTooManyNumbers() {
    PolylineParser pp(2);
    pp.parse("(1,2,3), (-1.2, 3.45, +6.79), (-1.2e-3, 3.45e19, -689.e-12)");
    const std::vector< Vec3 >& points = pp.getPoints();
    pp.print();
    assert(points.size() == 3);    
}

void testTooFewNumbers() {
    PolylineParser pp(3);
    pp.parse("(1,2), (-1.2, 3.45), (-1.2e-3, 3.45e19)");
    const std::vector< Vec3 >& points = pp.getPoints();
    pp.print();
    assert(points.size() == 3);
}


int main() {
    test3d();
    test2d();
    testTooManyNumbers();
    testTooFewNumbers();
}

