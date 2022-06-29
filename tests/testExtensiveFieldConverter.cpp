#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntExtensiveFieldConverter.h>
#include <mntVectorInterp.h>
#include <mntLogger.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test1Cell() {
    Grid_t* grd = NULL;
    mnt_grid_new(&grd);

    vtkIdType numCells = 1;

    std::vector<double> points({
        0., 0., 0.,
        1., 0., 0.,
        1., 1., 0.,
        0., 1., 0.,
    });
    mnt_grid_setPointsPtr(&grd, &points[0]);
    mnt_grid_build(&grd, MNT_NUM_VERTS_PER_QUAD, numCells);

    ExtensiveFieldConverter_t* efc = NULL;
    double aRadius = 1;
    int ier = mnt_extensivefieldconverter_new(&efc, aRadius);

    int degrees = 0;
    ier = mnt_extensivefieldconverter_setGrid(&efc, grd, degrees);
    assert(ier == 0);

    std::vector<double> data(4);

    double tol = 1.e-15;

    std::vector<double> uedge({1.0, 0.0, 3.0, 0.0});
    std::vector<double> vedge({0.0, 2.0, 0.0, 4.0});
    ier = mnt_extensivefieldconverter_getEdgeDataFromCellByCellVectors(&efc, &uedge[0], &vedge[0], &data[0]);
    std::cerr << "test1Cell: edge data = " << data[0] << ',' << data[1] << ',' << data[2] << ',' << data[3] << '\n';
    assert(ier == 0);
    assert(fabs(data[0] - (+1.0)) < tol);
    assert(fabs(data[1] - (+2.0)) < tol);
    assert(fabs(data[2] - (+3.0)) < tol);
    assert(fabs(data[3] - (+4.0)) < tol);

    std::vector<double> uface({0.0, 2.0, 0.0, 4.0});
    std::vector<double> vface({1.0, 0.0, 3.0, 0.0});
    ier = mnt_extensivefieldconverter_getFaceDataFromCellByCellVectors(&efc, &uface[0], &vface[0], &data[0]);
    std::cerr << "test1Cell: face data = " << data[0] << ',' << data[1] << ',' << data[2] << ',' << data[3] << '\n';
    assert(ier == 0);
    assert(fabs(data[0] - (+1.0)) < tol);
    assert(fabs(data[1] - (+2.0)) < tol);
    assert(fabs(data[2] - (+3.0)) < tol);
    assert(fabs(data[3] - (+4.0)) < tol);



    mnt_grid_del(&grd);
}


void testUgridData() {

    int ier;
    Grid_t* grd = NULL;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    /* a single cell
     3....>2....2
     :          :
     v          ^
     1          0
     :          :
     0....<3....1
    */

    std::size_t ncells = 1;
    std::size_t nedges = 4;
    std::size_t npoints = 4;
    std::vector<double> xyz{0.,0.,0.,
                            1.,0.,0.,
                            1.,1.,0.,
                            0.,1.,0.};
    std::vector<std::size_t> face2nodes{0, 1, 2, 3};
    std::vector<std::size_t> edge2nodes{1, 2,
                                        3, 0,
                                        3, 2,
                                        1, 0};
    ier = mnt_grid_loadFromUgrid2DData(&grd, ncells, nedges, npoints, &xyz[0], &face2nodes[0], &edge2nodes[0]);

    ExtensiveFieldConverter_t* efc = NULL;
    double aRadius = 1;
    ier = mnt_extensivefieldconverter_new(&efc, aRadius);

    int degrees = 0;
    ier = mnt_extensivefieldconverter_setGrid(&efc, grd, degrees);
    assert(ier == 0);

    std::vector<double> data(4);

    double tol = 1.e-15;

    std::vector<double> uedge({0.0, 0.0, 3.0, 1.0});
    std::vector<double> vedge({2.0, 4.0, 0.0, 0.0});
    ier = mnt_extensivefieldconverter_getEdgeDataFromUniqueEdgeVectors(&efc, &uedge[0], &vedge[0], &data[0]);
    std::cerr << "testUgridData: edge data = " << data[0] << ',' << data[1] << ',' << data[2] << ',' << data[3] << '\n';
    assert(ier == 0);
    assert(fabs(data[0] - (+1.0)) < tol);
    assert(fabs(data[1] - (+2.0)) < tol);
    assert(fabs(data[2] - (+3.0)) < tol);
    assert(fabs(data[3] - (+4.0)) < tol);

    std::vector<double> uface({2.0, 4.0, 0.0, 0.0});
    std::vector<double> vface({0.0, 0.0, 3.0, 1.0});
    ier = mnt_extensivefieldconverter_getFaceDataFromUniqueEdgeVectors(&efc, &uface[0], &vface[0], &data[0]);
    std::cerr << "testUgridData: face data = " << data[0] << ',' << data[1] << ',' << data[2] << ',' << data[3] << '\n';
    assert(ier == 0);
    assert(fabs(data[0] - (+1.0)) < tol);
    assert(fabs(data[1] - (+2.0)) < tol);
    assert(fabs(data[2] - (+3.0)) < tol);
    assert(fabs(data[3] - (+4.0)) < tol);

    // clean up
    ier = mnt_grid_del(&grd);
    ier = mnt_extensivefieldconverter_del(&efc);

}


void testUniformDegrees() {

    double aRadius = 1;

    int ier;
    Grid_t* grd = NULL;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    int fixLonAcrossDateline = 0; // lat-lon
    int averageLonAtPole = 0; // lat-lon
    int degrees = 1;
    mnt_grid_setFlags(&grd, fixLonAcrossDateline, averageLonAtPole, degrees);

    // filename$meshname
    mnt_grid_loadFromUgrid2DFile(&grd, "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc$latlon");

    std::size_t numCells = 0;
    mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(numCells > 0);

    std::size_t numEdges = 0;
    mnt_grid_getNumberOfEdges(&grd, &numEdges);
    assert(numEdges > 0);

    std::vector<double> uEdge(numEdges);
    std::vector<double> vEdge(numEdges);
    std::vector<double> dataEdge(numCells * MNT_NUM_VERTS_PER_QUAD);
    std::vector<double> dataFace(numCells * MNT_NUM_VERTS_PER_QUAD);

    std::size_t edgeId;
    int signEdge;
    std::vector<double> p0(3), p1(3);
    for (std::size_t cellId = 0; cellId < numCells; ++cellId) {
        for (int edgeIndex = 0; edgeIndex < MNT_NUM_VERTS_PER_QUAD; ++edgeIndex) {
            mnt_grid_getEdgeId(&grd, cellId, edgeIndex, &edgeId, &signEdge);
            uEdge[edgeId] = 1.0; // in m/s
            vEdge[edgeId] = 0.0; // in m/s
        }
    }

    ExtensiveFieldConverter_t* efc = NULL;
    ier = mnt_extensivefieldconverter_new(&efc, aRadius);

    ier = mnt_extensivefieldconverter_setGrid(&efc, grd, degrees);
    assert(ier == 0);

    ier = mnt_extensivefieldconverter_getEdgeDataFromUniqueEdgeVectors(&efc, &uEdge[0], &vEdge[0], &dataEdge[0]);
    assert(ier == 0);

    ier = mnt_extensivefieldconverter_getFaceDataFromUniqueEdgeVectors(&efc, &uEdge[0], &vEdge[0], &dataFace[0]);
    assert(ier == 0);

    for (auto cellId = 0; cellId < numCells; ++cellId) {
        for (auto edgeIndex = 0; edgeIndex < MNT_NUM_VERTS_PER_QUAD; ++edgeIndex) {
            auto k = cellId*MNT_NUM_VERTS_PER_QUAD + edgeIndex;
            std::cout << "testUniformDegrees: edge integrated data cellId = " << cellId << " edgeIndex = " << edgeIndex << " edge val = " << dataEdge[k] << " face val = " << dataFace[k] << '\n';
        }
    }

    VectorInterp_t* vi = NULL;
    mnt_vectorinterp_new(&vi);
    mnt_vectorinterp_setGrid(&vi, grd);
    int numCellsPerBucket = 128;
    double periodX = 360;
    int enableFolding = 0;
    mnt_vectorinterp_buildLocator(&vi, numCellsPerBucket, periodX, enableFolding);

    std::size_t numPoints = 1;
    std::vector<double> targetPoints({
        45., 0., 0.
    });
    double tol = 1.e-10;
    mnt_vectorinterp_findPoints(&vi, numPoints, &targetPoints[0], tol);

    // data is always defined cell by cell
    Vec3 targetVector;
    ier = mnt_vectorinterp_getEdgeVectorsFromCellByCellData(&vi, &dataEdge[0], &targetVector[0]);
    assert(ier == 0);
    for (auto it = 0; it < numPoints; ++it) {
        std::cout << "testUniformDegrees: target edge vector at point " << targetPoints[it*3+0] << ',' << targetPoints[it*3+1] << " is " << targetVector << '\n';
    }

    // clean up
    mnt_grid_del(&grd);
    mnt_vectorinterp_del(&vi);
    ier = mnt_extensivefieldconverter_del(&efc);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    testUniformDegrees();
    // test1Cell();
    // testUgridData();

    mnt_printLogMessages();

    return 0;
}
