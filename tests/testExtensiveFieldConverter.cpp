#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntExtensiveFieldConverter.h>
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


void testLFRic() {
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);
    // filename$meshname
    ier = mnt_grid_loadFromUgrid2DFile(&grd, "${CMAKE_SOURCE_DIR}/data/lfric_diag_wind.nc$Mesh2d");
    assert(ier == 0);
    std::size_t numBadCells = 0;
    ier = mnt_grid_check(&grd, &numBadCells);
    assert(numBadCells == 0);
    double* verts;
    ier = mnt_grid_getPointsPtr(&grd, &verts);
    assert(ier == 0);
    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    for (auto cellId = 0; cellId < numCells; ++cellId) {
        for (auto vIndex = 0; vIndex < 4; ++vIndex) {
            auto k = 3*(cellId*4 + vIndex);
            std::cout << "testLFRic: cell " << cellId << " node " << vIndex << " coords " << 
                         verts[k + 0] << ',' << verts[k + 1] << ',' << verts[k + 2] << '\n';
        }
    }
    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    // testLFRic();
    test1Cell();
    testUgridData();

    mnt_printLogMessages();

    return 0;
}
