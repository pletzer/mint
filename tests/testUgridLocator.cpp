#include <mntUgrid2D.h>
#include <iostream>
#include <cmath>
#undef NDEBUG // turn on asserts


void testPointInside() {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    /* 
    faceId 37 nodes
    180 -0 0 
    180 22.5 0 
    157.5 20.94102047224 0 
    157.5 -0 0 
    */

    const double tol = 1.e-10;
    const size_t faceId = 37;
    size_t nSegs = 10;
    const double pBegPtr[] = {180., 0., 0.};
    const double pEndPtr[] = {157.5, 20.94102047224, 0.};
    Vec3 pBeg(pBegPtr);
    Vec3 pEnd(pEndPtr);
    Vec3 du = pEnd - pBeg;
    du /= (double) nSegs;
    for (size_t i = 0; i < nSegs + 1; ++i) {
        Vec3 p = pBeg + (double) i * du;
        std::vector<Vec3> nodes = ug.getFacePointsRegularized(faceId);
        std::cerr << "checking if point " << p << " is inside of cell " << faceId << '\n';
        for (size_t j = 0; j < nodes.size(); ++j) {
            std::cout << nodes[j] << '\n';
        }
        assert(ug.containsPoint(faceId, p, tol));
    }

}

void testPointOutside() {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    /* 
    faceId 37 nodes
    180 -0 0 
    180 22.5 0 
    157.5 20.94102047224 0 
    157.5 -0 0 
    */

    const double tol = 1.e-10;
    const size_t faceId = 37;
    size_t nSegs = 10;
    const double pBegPtr[] = {190., 0., 0.};
    const double pEndPtr[] = {180.00001, 0., 0.};
    Vec3 pBeg(pBegPtr);
    Vec3 pEnd(pEndPtr);
    Vec3 du = pEnd - pBeg;
    du /= (double) nSegs;
    for (size_t i = 0; i < nSegs + 1; ++i) {
        Vec3 p = pBeg + (double) i * du;
        std::vector<Vec3> nodes = ug.getFacePointsRegularized(faceId);
        std::cerr << "checking if point " << p << " is outside of cell " << faceId << '\n';
        for (size_t j = 0; j < nodes.size(); ++j) {
            std::cout << nodes[j] << '\n';
        }
        assert(!ug.containsPoint(faceId, p, tol));
    }

}


void testPoint(const Vec3& p) {

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
    bool found = ug.findCell(p, tol, &faceId);
    if (found) {
        std::vector<Vec3> nodes = ug.getFacePointsRegularized(faceId);
        std::cout << "OK: point " << p << " is inside face " << faceId << " which has nodes:\n";
        for (const Vec3& node : nodes) {
            std::cout << node << '\n';
        }
        assert(ug.containsPoint(faceId, p, 1.e-8));
    }
    else {
        std::cout << "point " << p << " was not found\n";
    }
    assert(found);

    // find the parametric coordinates
    ug.setCellPoints(faceId);
    Vec3 pcoords;
    found = ug.getParamCoords(p, &pcoords[0]);
    assert(found);

    // check that the pcoords are correct
    Vec3 p2;
    ug.interpolate(pcoords, &p2[0]);

    Vec3 dp = p2 - p;
    double error = std::sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    std::cerr << "dist^2 error = " << error << '\n';
    assert(error < 1.e-6);
}

void testLine(const Vec3& p0, const Vec3& p1) {

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

    // check that we found all the cells by dividing the line in 1000 segments, checking that each point 
    // is inside of the cells we found
    size_t nSegments = 1000;
    Vec3 du = p1 - p0;
    du /= (double) nSegments;
    size_t cId;
    const double tol = 1.e-14;
    for (size_t iSegment = 0; iSegment < nSegments; ++iSegment) {
        Vec3 p = p0 + (double) iSegment * du;
        bool found = ug.findCell(p, tol, &cId);
        if (found) {
            if (faceIds.find(cId) == faceIds.end()) {
                std::cerr << "ERROR: unable to find point " << p << " belonging to face " << cId 
                          << " and along line " << p0 << " -> " << p1 << " among the above faces/cells!\n";
                assert(false);
            }
        }
    }
}

void testLineOutsideDomain() {

    std::string file = "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    // build the locator
    ug.buildLocator(1);

    const double p0Ptr[] = {1., 0., 0.};
    const double p1Ptr[] = {2., 0., 0.};
    Vec3 p0(p0Ptr);
    Vec3 p1(p1Ptr);
    size_t numFaces = ug.getNumberOfFaces();

    std::set<size_t> faceIds = ug.findCellsAlongLine(p0, p1);
    std::cout << "point " << p0 << " -> " << p1 << "overlaps with " << faceIds.size() << " cells:\n";
    for (const size_t& faceId : faceIds) {
        std::cout << faceId << ' ';
        assert(faceId < numFaces);
    }
    std::cout << '\n';

    // check that we found all the cells by dividing the line in 1000 segments, checking that each point 
    // is inside of the cells we found
    size_t nSegments = 1000;
    Vec3 du = p1 - p0;
    du /= (double) nSegments;
    size_t cId;
    const double tol = 1.e-14;
    for (size_t iSegment = 0; iSegment < nSegments; ++iSegment) {
        Vec3 p = p0 + (double) iSegment * du;
        bool found = ug.findCell(p, tol, &cId);
        if (found) {
            if (faceIds.find(cId) == faceIds.end()) {
                std::cerr << "ERROR: unable to find point " << p << " belonging to face " << cId 
                          << " and along line " << p0 << " -> " << p1 << " among the above faces/cells!\n";
                assert(false);
            }
        }
    }
}



int main() {

    testLineOutsideDomain();

    Vec3 p0(0.0);
    Vec3 p1(0.0);

    p0[0] =   0.; p0[1] = -67.;
    p1[0] = 360.; p1[1] =  67.;
    testLine(p0, p1);

    // going the other way
    p0[0] = 360.; p0[1] =  67.;
    p1[0] =   0.; p1[1] = -67.;
    testLine(p0, p1);

    testPointInside();
    testPointOutside();

    p0[0] = 180.; p0[1] = 0.;
    testPoint(p0);

    p0[0] = 90.; p0[1] = 90.;
    testPoint(p0);

    p0[0] = 90.; p0[1] = -90.;
    testPoint(p0);

    p0[0] = 0.; p0[1] = 67.;
    testPoint(p0);


    return 0;
}   
