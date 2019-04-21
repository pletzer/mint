#include <mntUgrid2D.h>
#include <iostream>
#undef NDEBUG // turn on asserts


void testPoint(const Vector<double>& p) {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    ug.dumpGridVtk("ug_cs_4.vtk");

    // build the locator
    ug.buildLocator(1);

    // find face
    const double tol = 1.e-10;
    size_t faceId;
    bool found = ug.findCell(p, tol, faceId);
    if (found) {
        std::vector< Vector<double> > nodes = ug.getFacePointsRegularized(faceId);
        std::cout << "OK: point " << p << " is inside face " << faceId << " which has nodes:\n";
        for (const Vector<double>& node : nodes) {
            std::cout << node << '\n';
        }
    }
    else {
        std::cout << "point " << p << " was not found\n";
    }
    assert(found);

}

void testLine(const Vector<double>& p0, const Vector<double>& p1) {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    // build the locator
    ug.buildLocator(1);

    std::set<size_t> faceIds = ug.findCellsAlongLine(p0, p1);
    std::cout << "point " << p0 << " -> " << p1 << "overlaps with " << faceIds.size() << " cells:\n";
    for (const size_t& faceId : faceIds) {
        std::cout << faceId << ' ';
    }
    std::cout << '\n';

}


int main() {

    Vector<double> p0(3, 0.0);

    p0[0] = 180.; p0[1] = 0.;
    testPoint(p0);

    p0[0] = 90.; p0[1] = 90.;
    testPoint(p0);

    p0[0] = 90.; p0[1] = -90.;
    testPoint(p0);

    p0[0] = 0.; p0[1] = 67.;
    testPoint(p0);

    Vector<double> p1(3, 0.0);
    p0[0] =   0.; p0[1] = -67.;
    p1[0] = 360.; p1[1] =  67.;
    testLine(p0, p1);


    return 0;
}   
