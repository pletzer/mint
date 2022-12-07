#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVectorInterp.h>
#include <mntMatMxN.h>
#include "saveEdgeVectors.h"

#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>

/**
 * Create uniform VTK grid objects
 * @param nx number of x cells
 * @param ny number of y cells
 * @param points array of size (4*ny*nx) to be filled in
 */
void createUniformGrid(int nx, int ny, double points[]) {

    double dx = 360. / (double) nx;
    double dy = 180. / (double) ny;

    int k = 0;
    for (int i = 0; i < nx; ++i) {
        double x0 = i*dx;
        double x1 = x0 + dx;
        for (int j = 0; j < ny; ++j) {
            double y0 = -90.0 + j*dy;
            double y1 = y0 + dy;

            points[4*k*3 +  0] = x0;
            points[4*k*3 +  1] = y0;
            points[4*k*3 +  2] = 0.0;

            points[4*k*3 +  3] = x1;
            points[4*k*3 +  4] = y0;
            points[4*k*3 +  5] = 0.0;

            points[4*k*3 +  6] = x1;
            points[4*k*3 +  7] = y1;
            points[4*k*3 +  8] = 0.0;

            points[4*k*3 +  9] = x0;
            points[4*k*3 + 10] = y1;
            points[4*k*3 + 11] = 0.0;

            k++;
        }
    }
}

void testSimple() {

    int ier;
    vtkIdType numCells = 1;
    std::vector<double> points{0., 0., 0.,
                               1., 0., 0.,
                               1., 1., 0.,
                               0., 1., 0.};
    Grid_t* mgrid = NULL;
    ier = mnt_grid_new(&mgrid);
    assert(ier == 0);
    int fixLonAcrossDateline = 0;
    int averageLonAtPole = 0;
    int degrees = 0;
    ier = mnt_grid_setFlags(&mgrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);
    ier = mnt_grid_setPointsPtr(&mgrid, &points[0]);
    assert(ier == 0);
    ier = mnt_grid_build(&mgrid, 4, numCells);
    assert(ier == 0);
    ier = mnt_grid_print(&mgrid);
    assert(ier == 0);

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    vtkUnstructuredGrid* ugrid;
    ier = mnt_grid_get(&mgrid, &ugrid);
    assert(ier == 0);
    cloc->SetDataSet(ugrid);
    cloc->SetNumberOfCellsPerBucket(10);
    cloc->BuildLocator();

    VectorInterp_t* vp = NULL;
    ier = mnt_vectorinterp_new(&vp);
    assert(ier == 0);
    ier = mnt_vectorinterp_setGrid(&vp, mgrid);
    assert(ier == 0);
    ier = mnt_vectorinterp_setLocator(&vp, cloc);
    assert(ier == 0);
    double tol2 = 1.e-10;

    // set the data
    std::vector<double> data{0., 1., 2., 3.};
    Vec3 target, resVec;

    target[0] = 0.;
    target[1] = 0.;
    target[2] = 0.;
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA,
                                                            &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0.0) < 1.e-10 && fabs(resVec[1] - 3.0) < 1.e-10);

    target[0] = 1.;
    target[1] = 0.;
    target[2] = 0.;
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA,
                                                            &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0.0) < 1.e-10 && fabs(resVec[1] - 1.0) < 1.e-10);

    target[0] = 1.;
    target[1] = 1.;
    target[2] = 0.;
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA,
                                                            &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 2.0) < 1.e-10 && fabs(resVec[1] - 1.0) < 1.e-10);

    target[0] = 0.;
    target[1] = 1.;
    target[2] = 0.;
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA,
                                                            &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 2.0) < 1.e-10 && fabs(resVec[1] - 3.0) < 1.e-10);

    // destroy
    ier = mnt_vectorinterp_del(&vp);
    assert(ier == 0);
    cloc->Delete();
    ier = mnt_grid_del(&mgrid);
    assert(ier == 0);
}

void testRotated() {

    int ier;
    vtkIdType numCells = 1;
    double a = 1.0;
    double b = 1.0;
    double angleDeg = 30.0;
    double cosAngle = cos(angleDeg*M_PI/180.);
    double sinAngle = sin(angleDeg*M_PI/180.);
    Mat3x3 rotation;
    rotation(0, 0) = cosAngle; rotation(0, 1) = -sinAngle; rotation(0, 2) = 0;
    rotation(1, 0) = sinAngle; rotation(1, 1) =  cosAngle; rotation(1, 2) = 0;
    rotation(2, 0) =        0; rotation(2, 1) =         0; rotation(2, 2) = 1;
    std::vector<double> pointsOri{0, 0, 0,
                                  a, 0, 0,
                                  a, b, 0,
                                  0, b, 0};
    std::vector<double> points(4*3);
    for (std::size_t i = 0; i < 4; ++i) {
        Vec3 p0;
        for (std::size_t j = 0; j < 3; ++j) {
            p0[j] = pointsOri[i*3 + j];
        }
        // rotate
        Vec3 p1 = dot(rotation, p0);
        for (std::size_t j = 0; j < 3; ++j) {
            points[i*3 + j] = p1[j];
        }
    }

    Grid_t* mgrid = NULL;
    ier = mnt_grid_new(&mgrid);
    assert(ier == 0);
    int fixLonAcrossDateline = 0;
    int averageLonAtPole = 0;
    int degrees = 0;
    ier = mnt_grid_setFlags(&mgrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);
    ier = mnt_grid_setPointsPtr(&mgrid, &points[0]);
    assert(ier == 0);
    ier = mnt_grid_build(&mgrid, 4, numCells);
    assert(ier == 0);
    ier = mnt_grid_print(&mgrid);
    assert(ier == 0);

    VectorInterp_t* vp = NULL;
    ier = mnt_vectorinterp_new(&vp);
    assert(ier == 0);
    ier = mnt_vectorinterp_setGrid(&vp, mgrid);
    assert(ier == 0);
    int numCellsPerBucket = 10;
    double periodX = 0.;
    int enableFolding = 0;
    ier = mnt_vectorinterp_buildLocator(&vp, numCellsPerBucket, periodX, enableFolding);
    double tol2 = 1.e-10;

    // set some data
    std::vector<double> data{1, 0, 0, 0};
    Vec3 target, resVec;

    // basis vector at south edge
    target[0] = points[0];
    target[1] = points[1];
    target[2] = points[2];
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0.8660254037844) < 1.e-12 && fabs(resVec[1] - 0.5) < 1.e-12);
    ier = mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " face vector is " << resVec << '\n';
    // note negative sign because area is pointing down
    assert(fabs(resVec[1] - (-0.8660254037844)) < 1.e-12 && fabs(resVec[0] - 0.5) < 1.e-12);

    target[0] = points[3];
    target[1] = points[4];
    target[2] = points[5];
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0.8660254037844) < 1.e-12 && fabs(resVec[1] - 0.5) < 1.e-12);
    ier = mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " face vector is " << resVec << '\n';
    // negative because area is pointing down
    assert(fabs(resVec[1] - (-0.8660254037844)) < 1.e-12 && fabs(resVec[0] - 0.5) < 1.e-12);

    target[0] = points[6];
    target[1] = points[7];
    target[2] = points[8];
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0) < 1.e-12 && fabs(resVec[1] - 0) < 1.e-12);

    target[0] = points[9];
    target[1] = points[10];
    target[2] = points[11];
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0) < 1.e-12 && fabs(resVec[1] - 0) < 1.e-12);

    // basis vector at east edge
    data[0] = 0.0; data[1] = 1.0; data[2] = 0.0; data[3] = 0.0;
    target[0] = 0.5*(points[3] + points[6]);
    target[1] = 0.5*(points[4] + points[7]);
    target[2] = 0.5*(points[5] + points[8]);
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " east edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] + 0.5) < 1.e-12 && fabs(resVec[1] - 0.8660254037844) < 1.e-12);
    ier = mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " east face vector is " << resVec << '\n';
    assert(fabs(resVec[1] - 0.5) < 1.e-12 && fabs(resVec[0] - 0.8660254037844) < 1.e-12);

    // basis vector at west edge
    data[0] = 0.0; data[1] = 0.0; data[2] = 0.0; data[3] = 1.0;
    target[0] = 0.5*(points[ 9] + points[0]);
    target[1] = 0.5*(points[10] + points[1]);
    target[2] = 0.5*(points[11] + points[2]);
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " west edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] + 0.5) < 1.e-12 && fabs(resVec[1] - 0.8660254037844) < 1.e-12);
    ier = mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " west face vector is " << resVec << '\n';
    assert(fabs(resVec[1] - 0.5) < 1.e-12 && fabs(resVec[0] - 0.8660254037844) < 1.e-12);

    // basis vector at north edge
    data[0] = 0.0; data[1] = 0.0; data[2] = 1.0; data[3] = 0.0;
    target[0] = 0.5*(points[ 6] + points[ 9]);
    target[1] = 0.5*(points[ 7] + points[10]);
    target[2] = 0.5*(points[ 8] + points[11]);
    ier = mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    assert(ier == 0);
    ier = mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " north edge vector is " << resVec << '\n';
    assert(fabs(resVec[0] - 0.8660254037844) < 1.e-12 && fabs(resVec[1] - 0.5) < 1.e-12);
    ier = mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    assert(ier == 0);
    std::cout << "at point " << target << " north face vector is " << resVec << '\n';
    // note area is pointing down
    assert(fabs(resVec[1] - (-0.8660254037844)) < 1.e-12 && fabs(resVec[0] - 0.5) < 1.e-12);

    // destroy
    ier = mnt_vectorinterp_del(&vp);
    assert(ier == 0);
    ier = mnt_grid_del(&mgrid);
    assert(ier == 0);
}


void testUniformGrid(int nx, int ny) {

    std::vector<double> points(4*ny*nx*3);
    createUniformGrid(nx, ny, &points[0]);

    Grid_t* mgrid = NULL;
    mnt_grid_new(&mgrid);
    int fixLonAcrossDateline = 0;
    int averageLonAtPole = 0;
    int degrees = 1;
    mnt_grid_setFlags(&mgrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    mnt_grid_setPointsPtr(&mgrid, &points[0]);
    mnt_grid_build(&mgrid, 4, ny*nx);

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    vtkUnstructuredGrid* ugrid;
    mnt_grid_get(&mgrid, &ugrid);
    cloc->SetDataSet(ugrid);
    cloc->SetNumberOfCellsPerBucket(10);
    cloc->BuildLocator();
    cloc->setPeriodicityLengthX(360.0);

    // set some data
    std::vector<double> data(4*ny*nx);
    for (std::size_t i = 0; i < 4*ny*nx; ++i) {
        data[i] = (double) i;
    }

    Vec3 target, resVec;
    // inside the first face
    target[0] = 0.00001;
    target[1] = -89.99999;
    target[2] = 0.;

    VectorInterp_t* vp = NULL;
    mnt_vectorinterp_new(&vp);
    mnt_vectorinterp_setGrid(&vp, mgrid);
    mnt_vectorinterp_setLocator(&vp, cloc);
    double tol2 = 1.e-10;
    mnt_vectorinterp_findPoints(&vp, 1, &target[0], tol2);
    mnt_vectorinterp_getEdgeVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    std::cout << "at point " << target << " edge vector is " << resVec << '\n';
    mnt_vectorinterp_getFaceVectors(&vp, &data[0], MNT_CELL_BY_CELL_DATA, &resVec[0]);
    std::cout << "at point " << target << " face vector is " << resVec << '\n';

    // destroy
    mnt_vectorinterp_del(&vp);
    cloc->Delete();
    mnt_grid_del(&mgrid);
}

void testCubedSphereFaceVectorsOnEdge() {

    int ier;

    // read the grid
    Grid_t* grd = NULL;
    mnt_grid_new(&grd);
    mnt_grid_setFlags(&grd, 1, 1, 1); // cubed-sphere, degrees
    ier = mnt_grid_loadFromUgrid2DFile(&grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc$physics");
    assert(ier == 0);

    std::size_t numEdges;
    mnt_grid_getNumberOfEdges(&grd, &numEdges);
    std::size_t numCells;
    mnt_grid_getNumberOfCells(&grd, &numCells);
    std::cout << "testCubedSphereFaceVectorsOnEdge: num cells: " << numCells << " num edges: " << numEdges << '\n';

    VectorInterp_t* vp = NULL;
    ier = mnt_vectorinterp_new(&vp);
    assert(ier == 0);
    ier = mnt_vectorinterp_setGrid(&vp, grd);
    assert(ier == 0);
    // mnt_vectorinterp_buildLocator(&vp, 256, 360., 0);

    std::vector<double> faceData(numEdges);
    std::vector<double> uedge(numEdges);
    std::vector<double> vedge(numEdges);
    std::vector<double> uedgeExact(numEdges);
    std::vector<double> vedgeExact(numEdges);
    std::size_t edgeId;
    int edgeSign;
    Vec3 p0, p1;
    std::vector<Vec3> p0s(numEdges);
    std::vector<Vec3> p1s(numEdges);

    for (std::size_t icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            ier = mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &edgeSign);
            assert(ier == 0);
            ier = mnt_grid_getPoints(&grd, icell, ie, &p0[0], &p1[0]);
            assert(ier == 0);
            double lam0 = p0[0] * M_PI/180;
            double lam1 = p1[0] * M_PI/180;
            double the0 = p0[1] * M_PI/180;
            double the1 = p1[1] * M_PI/180;
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            double lamMid = 0.5*(lam0 + lam1);
            double theMid = 0.5*(the0 + the1);

            faceData[edgeId] = edgeSign * (s1 - s0);
            uedgeExact[edgeId] = - sin(theMid) * cos(lamMid);
            vedgeExact[edgeId] = sin(lamMid);

            p0s[edgeId] = p0;
            p1s[edgeId] = p1;
        }
    }

    ier = mnt_vectorinterp_getVectorsOnEdges(&vp,
        &faceData[0], MNT_UNIQUE_EDGE_DATA, &uedge[0], &vedge[0], MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    // check
    double error = 0;
    for (auto i = 0; i < numEdges; ++i) {
        error += fabs(uedgeExact[i] - uedge[i]) + fabs(vedgeExact[i] - vedge[i]);
        // std::cerr << "testCubedSphereFaceVectorsOnEdge: edge=" << i << " points " << p0s[edgeId] << ';' << p1s[edgeId]
        //           << " u=" << uedge[i] << " (exact " << uedgeExact[i] 
        //           << ") v=" << vedge[i] << " (exact " << vedgeExact[i] << ")\n";
    }
    error /= numEdges;
    std::cerr << "testCubedSphereFaceVectorsOnEdge error: " << error << '\n';

    mnt_grid_dump(&grd, "testCubedSphereFaceVectorsOnEdge_grid.vtk");
    saveEdgeVectors(grd, uedge, vedge, "testCubedSphereFaceVectorsOnEdge_vectors.vtk");
    saveEdgeVectors(grd, uedgeExact, vedgeExact, "testCubedSphereFaceVectorsOnEdge_vectorsExact.vtk");
    assert(error < 0.005);


    // destroy
    mnt_vectorinterp_del(&vp);
    mnt_grid_del(&grd);
}

void testSqueezed() {
    int ier;

    // create grid
    Grid_t* grid = NULL;
    mnt_grid_new(&grid);
    mnt_grid_setFlags(&grid, 0, 0, 1);
    std::size_t ncells = 1, nedges = 4, npoints = 4;
    const double xyz[] = {
        135., 90., 0.,
        90., 84.375, 0.,
        135., 82.07041073676, 0.,
        180., 84.375, 0. 
    };
    const std::size_t face2nodes[] = {
        0, 1, 2, 3
    };
    const std::size_t edge2nodes[] = {
        0, 1,
        1, 2,
        3, 2,
        0, 3
    };
    ier = mnt_grid_loadFromUgrid2DData(&grid, ncells, nedges, npoints, 
                                 xyz, face2nodes, edge2nodes);
    assert(ier == 0);
    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grid, &numCells);
    assert(ier == 0);
    std::size_t numEdges;
    ier = mnt_grid_getNumberOfEdges(&grid, &numEdges);

    Vec3 p0, p1, pMid;

    std::vector<double> uCellByCellExact(numCells * MNT_NUM_EDGES_PER_QUAD);
    std::vector<double> vCellByCellExact(numCells * MNT_NUM_EDGES_PER_QUAD);
    for (auto icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            ier = mnt_grid_getPoints(&grid, icell, ie, &p0[0], &p1[0]);
            assert(ier == 0);
            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + ie;
            pMid = 0.5*(p0 + p1);
            uCellByCellExact[k] = -sin(pMid[1] * M_PI/180.) * cos(pMid[0] * M_PI/180.);
            vCellByCellExact[k] = sin(pMid[0] * M_PI/180.);
        }
    }

    std::vector<double> fluxes({0., -0.09755, 0.000467, -0.098017});
    VectorInterp_t* vp = NULL;
    ier = mnt_vectorinterp_new(&vp);
    assert(ier == 0);
    ier = mnt_vectorinterp_setGrid(&vp, grid);
    assert(ier == 0);
    std::vector<double> uEdges(numEdges);
    std::vector<double> vEdges(numEdges);
    ier = mnt_vectorinterp_getVectorsOnEdges(&vp, &fluxes[0],
                           MNT_CELL_BY_CELL_DATA, &uEdges[0], &vEdges[0], MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    // check
    for (auto i = 0; i < 4; ++i) {
        std::cout << "testSqueezed: u=" << uEdges[i] << " exact=" << uCellByCellExact[i] << " diff: " << uCellByCellExact[i] - uEdges[i] << '\n';
        assert(std::fabs(uEdges[i] - uCellByCellExact[i]) < 0.012);
        std::cout << "testSqueezed: v=" << vEdges[i] << " exact=" << vCellByCellExact[i] << " diff: " << vCellByCellExact[i] - vEdges[i] << '\n';
        assert(std::fabs(vEdges[i] - vCellByCellExact[i]) < 0.012);
    }


    // clean up
    ier = mnt_vectorinterp_del(&vp);
    assert(ier == 0);
    ier = mnt_grid_del(&grid);
    assert(ier == 0);

}

int main(int argc, char** argv) {

    testSqueezed();
    testCubedSphereFaceVectorsOnEdge();
    testSimple();
    testRotated();
    testUniformGrid(8, 4);

    return 0;
}
