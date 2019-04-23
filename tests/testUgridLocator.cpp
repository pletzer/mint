#include <mntUgrid2D.h>
#include <iostream>
#include <cmath>
#undef NDEBUG // turn on asserts

void testLineGridIntersections(const Vector<double>& pBeg,
                               const Vector<double>& pEnd) {

    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ug;

    int ier = ug.load(file, "physics");
    assert(ier == 0);

    ug.buildLocator(10);

    std::cout << "testLineGridIntersections: " << pBeg << " -> " << pEnd << '\n';

    std::vector< std::pair<size_t, std::vector<double> > > 
                      cIdlambdas = ug.findIntersectionsWithLine(pBeg, pEnd);

    Vector<double> u = pEnd - pBeg;

    // check that all the lambda intervals add to one
    double totLambda = 0;
    printf("  cell      lambdaBeg       lambdaEnd\n");
    for (const std::pair< size_t, std::vector<double> >& kv : cIdlambdas) {
        size_t cellId = kv.first;
        double lambdaBeg = kv.second[0];
        double lambdaEnd = kv.second[1];
        printf("%6ld     %10.6lf      %10.6lf\n", cellId, lambdaBeg, lambdaEnd);
        totLambda += lambdaEnd - lambdaBeg;
        // make sure the start and end points are inside the cell
        Vector<double> p0 = pBeg + lambdaBeg*u;
        assert(ug.containsPoint(cellId, p0, 1.e-10));
        Vector<double> p1 = pBeg + lambdaEnd*u;
        assert(ug.containsPoint(cellId, p1, 1.e-10));    
    }
    std::cout << "total integrated lambda: " << totLambda << '\n';
    assert(std::abs(totLambda - 1.) < 1.e-10);
}

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
    Vector<double> pBeg{180., 0., 0.};
    Vector<double> pEnd{157.5, 20.94102047224, 0.};
    Vector<double> du = pEnd - pBeg;
    du /= (double) nSegs;
    for (size_t i = 0; i < nSegs + 1; ++i) {
        Vector<double> p = pBeg + (double) i * du;
        std::vector< Vector<double> > nodes = ug.getFacePointsRegularized(faceId);
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
    Vector<double> pBeg{190., 0., 0.};
    Vector<double> pEnd{180.00001, 0., 0.};
    Vector<double> du = pEnd - pBeg;
    du /= (double) nSegs;
    for (size_t i = 0; i < nSegs + 1; ++i) {
        Vector<double> p = pBeg + (double) i * du;
        std::vector< Vector<double> > nodes = ug.getFacePointsRegularized(faceId);
        std::cerr << "checking if point " << p << " is outside of cell " << faceId << '\n';
        for (size_t j = 0; j < nodes.size(); ++j) {
            std::cout << nodes[j] << '\n';
        }
        assert(!ug.containsPoint(faceId, p, tol));
    }

}


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
    bool found = ug.findCell(p, tol, &faceId);
    if (found) {
        std::vector< Vector<double> > nodes = ug.getFacePointsRegularized(faceId);
        std::cout << "OK: point " << p << " is inside face " << faceId << " which has nodes:\n";
        for (const Vector<double>& node : nodes) {
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
    Vector<double> pcoords(3);
    found = ug.getParamCoords(p, &pcoords[0]);
    assert(found);

    // check that the pcoords are correct
    Vector<double> p2(3);
    ug.interpolate(pcoords, &p2[0]);

    Vector<double> dp = p2 - p;
    double error = std::sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
    std::cerr << "dist^2 error = " << error << '\n';
    assert(error < 1.e-6);
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

    // check that we found all the cells by dividing the line in 1000 segments, checking that each point 
    // is inside of the cells we found
    size_t nSegments = 1000;
    Vector<double> du = p1 - p0;
    du /= (double) nSegments;
    size_t cId;
    const double tol = 1.e-14;
    for (size_t iSegment = 0; iSegment < nSegments; ++iSegment) {
        Vector<double> p = p0 + (double) iSegment * du;
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

    // because of the coarse cubed-sphere dicretization some paths may 
    // fall out of the domain even when longitude is in [0, 360] and 
    // latitude is in [-90, 90]. At some places 67 deg is the max latitude.

    Vector<double> p0(3, 0.0);
    Vector<double> p1(3, 0.0);

    p0[0] = 360.; p0[1] =  67.;
    p1[0] =   0.; p1[1] = -67.;
    // this is currently failing!!!
    //testLineGridIntersections(p0, p1);

    p0[0] =  90.; p0[1] = -90.;
    p1[0] =  90.; p1[1] = -90.;
    testLineGridIntersections(p0, p1);    

    p0[0] =  90.; p0[1] = -90.;
    p1[0] = 270.; p1[1] =  90.;
    testLineGridIntersections(p0, p1);    

    p0[0] =  90.; p0[1] = -90.;
    p1[0] =  90.; p1[1] =  90.;
    testLineGridIntersections(p0, p1);    

    p0[0] =  12.3; p0[1] = -67.;
    p1[0] =  12.3; p1[1] =  67.;
    testLineGridIntersections(p0, p1);    

    p0[0] =   0.; p0[1] = -67.;
    p1[0] = 360.; p1[1] =  67.;
    testLineGridIntersections(p0, p1);

    p0[0] =   0.; p0[1] = 0.;
    p1[0] = 360.; p1[1] = 0.;
    testLineGridIntersections(p0, p1);    

    p0[0] =   0.; p0[1] = -67.;
    testPoint(p0);
    p1[0] =  10.; p1[1] = -50.;
    testPoint(p1);
    testLineGridIntersections(p0, p1);

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

    p0[0] =   0.; p0[1] = -67.;
    p1[0] = 360.; p1[1] =  67.;
    testLine(p0, p1);


    return 0;
}   
