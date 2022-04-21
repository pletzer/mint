#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntPolylineIntegral.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

void test2Triangles() {
    
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    /* a single cell
     3....>2....2
     : \        :
     v  \       ^
     1   \      0
     :    ^     :
     :     4    :
     :      \   :
     :       \  :
     :        \ :
     :         \:
     0....<3....1
    */

    std::size_t ncells = 2;
    std::size_t nedges = 6;
    std::size_t npoints = 4;
    std::vector<double> xyz{0.,0.,0.,
                            1.,0.,0.,
                            1.,1.,0.,
                            0.,1.,0.};
    std::vector<std::size_t> face2nodes{0, 1, 3, 3,  // degenerate node
                                        1, 2, 3, 3}; // degenerate node
    std::vector<std::size_t> edge2nodes{1, 2,
                                        3, 0,
                                        3, 2,
                                        1, 0,
                                        1, 3,
                                        3, 3}; // degenerate node
    ier = mnt_grid_loadFromUgrid2DData(&grd, ncells, nedges, npoints, &xyz[0], &face2nodes[0], &edge2nodes[0]);
    assert(ier == 0);

    PolylineIntegral_t* pli = NULL;
    ier = mnt_polylineintegral_new(&pli);
    assert(ier == 0);

    ier = mnt_polylineintegral_setGrid(&pli, grd);
    assert(ier == 0);

    int numCellsPerBucket = 128;
    double periodX = 0;
    int enableFolding = 0;
    ier = mnt_polylineintegral_buildLocator(&pli, numCellsPerBucket, periodX, enableFolding);
    assert(ier == 0);

    int counterclock = 0;
    ier = mnt_polylineintegral_computeWeights(&pli, npoints, (const double*) &xyz[0], 
                                              counterclock);
    assert(ier == 0);

    std::vector<double> data{0, 1, 2, 3, 4, 5};
    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0],
                                           MNT_CELL_BY_CELL_DATA, &totalFlux);
    assert(ier == 0);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    test2Triangles();

    return 0;
}
