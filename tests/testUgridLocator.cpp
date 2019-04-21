#include <mntUgrid2D.h>
#include <iostream>
#undef NDEBUG // turn on asserts


void testPoint(const Vector<double>& p) {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    // build the locator
    ug.buildLocator(1);

    // find face
    const double tol = 1.e-10;
    size_t faceId;
    bool found = ug.findCell(p, tol, faceId);
    if (found) {
        std::vector< Vector<double> > nodes = ug.getFacePointsRegularized(faceId);
        std::cout << "point " << p << " is inside face " << faceId << " which has nodes:\n";
        for (const Vector<double>& node : nodes) {
            std::cout << node << '\n';
        }
    }
    else {
        std::cout << "point " << p << " was not found\n";
    }
    assert(found);

}


int main() {

    Vector<double> p(3, 0.0);

    p[0] = 180.; p[1] = 0.;
    testPoint(p);

    //p[0] = 180.; p[1] = 89.;
    //testPoint(p);

    return 0;
}   
