#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntPolylineIntegral.h>
#include <mntLogger.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

double potentialFunc(const double* p) {
    return p[0] + p[1];
}


void test1Triangle() {
    
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    /* 1 triangle
     3,2
     : \
     v  \
     3   \
     :    ^
     :     1
     :      \
     :       \ 
     :        \ 
     :         \
     0....0>....1
    */

    std::size_t ncells = 1;
    std::size_t nedges = 4;
    std::size_t npoints = 4;
    std::vector<double> xyz{0.,0.,0.,
                            1.,0.,0.,
                            1.,1.,0., // 0.,1.,0.,
                            0.,1.,0.}; // degenerate
    std::vector<std::size_t> face2nodes{0, 1, 2, 3};
    std::vector<std::size_t> edge2nodes{0, 1,
                                        1, 2,
                                        3, 2,
                                        0, 3};
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

    std::vector<double> targetPoints{0., 0., 0.,
                                     1., 0., 0.};
    std::size_t numTargetPoints = targetPoints.size() / 3;

    int counterclock = 0;
    ier = mnt_polylineintegral_computeWeights(&pli, numTargetPoints, (const double*) &targetPoints[0], 
                                              counterclock);
    mnt_printLogMessages();
    assert(ier == 0);

    // set the data on each edge
    std::vector<double> data{10, 11, 0, 13}; // set the flux to zero for the degenerate edge

    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0],
                                           MNT_UNIQUE_EDGE_DATA, &totalFlux);
    assert(ier == 0);

    // check
    double exactFlux = 10;
    double error = totalFlux - exactFlux;
    std::cout << "total flux: " << totalFlux << " exact flux: " << exactFlux << " error: " << error << '\n';
    // mnt_printLogMessages();
    assert(std::abs(error) < 1.e-10);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}


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

    // set the data on each edge
    std::vector<double> data(nedges);
    for (std::size_t iedge = 0; iedge < nedges; ++iedge) {

        std::size_t i0 = edge2nodes[iedge*MNT_NUM_VERTS_PER_EDGE + 0];
        const double* p0 = &xyz[i0*MNT_NUM_VERTS_PER_QUAD];
        double pot0 = potentialFunc(p0);

        std::size_t i1 = edge2nodes[iedge*MNT_NUM_VERTS_PER_EDGE + 1];
        const double* p1 = &xyz[i1*MNT_NUM_VERTS_PER_QUAD];
        double pot1 = potentialFunc(p1);

        data[iedge] = pot1 - pot0;
    }

    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0],
                                           MNT_UNIQUE_EDGE_DATA, &totalFlux);
    assert(ier == 0);

    // check
    double exactFlux = potentialFunc(&xyz[(npoints - 1)*3]) - potentialFunc(&xyz[0]);
    double error = totalFlux - exactFlux;
    std::cout << "total flux: " << totalFlux << " exact flux: " << exactFlux << " error: " << error << '\n';
    mnt_printLogMessages();
    assert(std::abs(error) < 1.e-10);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}


int main(int argc, char** argv) {

    test1Triangle();
    // test2Triangles();

    return 0;
}
