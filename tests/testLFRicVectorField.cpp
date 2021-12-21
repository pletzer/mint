#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntVectorInterp.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>


void testZonal() {
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);
    // one or more "$" to discriminate file and mesh names
    ier = mnt_grid_loadFromUgrid2D(&grd, "${CMAKE_SOURCE_DIR}/data/lfric_diag_wind.nc$Mesh2d");
    assert(ier == 0);
    std::size_t numBadCells = 0;
    ier = mnt_grid_check(&grd, &numBadCells);
    assert(numBadCells == 0);
 

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
    mnt_printLogMessages();
}

int main(int argc, char** argv) {

    testZonal();

    return 0;
}
