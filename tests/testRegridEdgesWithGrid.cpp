#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntRegridEdges.h>
#include <mntLogger.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test1() {

    int ier = 0, degrees = 1;
    int fixLonAcrossDateline, averageLonAtPole;

    // source grid
    Grid_t* srcGrid = NULL;
    ier = mnt_grid_new(&srcGrid);
    assert(ier == 0);
    // lon-lat grid does not require pole averaging
    fixLonAcrossDateline = 0;
    averageLonAtPole = 0;
    ier = mnt_grid_setFlags(&srcGrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);
    // read the grid
    ier = mnt_grid_loadFromUgrid2D(&srcGrid, "${CMAKE_SOURCE_DIR}/data/lonlatzt_100x50x3x2.nc$mesh2d");
    assert(ier == 0);

    // destination grid
    Grid_t* dstGrid = NULL;
    ier = mnt_grid_new(&dstGrid);
    assert(ier == 0);
    // cubed-sphere requires pole averaging
    fixLonAcrossDateline = 1;
    averageLonAtPole = 1;
    ier = mnt_grid_setFlags(&dstGrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);
    // read the grid
    ier = mnt_grid_loadFromUgrid2D(&dstGrid, "${CMAKE_SOURCE_DIR}/data/lfric_diag_wind.nc$Mesh2d");
    assert(ier == 0);

    // regridder
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, srcGrid);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dstGrid);
    assert(ier == 0);
    int numCellsPerBucket = 100;
    double periodX = 360.;
    int enableFolding = 0;
    ier = mnt_regridedges_buildLocator(&rgd, numCellsPerBucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dstGrid);
    mnt_grid_del(&srcGrid);
}

int main(int argc, char** argv) {

    test1();

    return 0;
}
