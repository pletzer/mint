#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntPolylineIntegral.h>
#include <mntLogger.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

#define EPS 1.e-5

double potentialFunc(const double* p) {
    return p[0] + p[1];
}


void test1Triangle(const std::vector<double>& targetLine, double exactFlux) {
    
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    /* 1 triangle
     3,2
     : \
     ^  \
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
                            EPS,1.,0.,
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

    std::size_t numTargetPoints = targetLine.size() / 3;

    int counterclock = 0;
    ier = mnt_polylineintegral_computeWeights(&pli, numTargetPoints, (const double*) &targetLine[0], 
                                              counterclock);
    assert(ier == 0);

    // set the data on each edge
    std::vector<double> data{10, 11, 0, 13}; // set the flux to zero for the degenerate edge

    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0],
                                           MNT_UNIQUE_EDGE_DATA, &totalFlux);
    assert(ier == 0);

    // check
    double error = totalFlux - exactFlux;
    std::cout << "total flux: " << totalFlux << " exact flux: " << exactFlux << " error: " << error << '\n';
    assert(std::abs(error) < 1.e-10);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    {
        std::vector<double> targetLine{0., 0., 0.,
                                       1., 0., 0.};
        double exactFlux = 10;
        test1Triangle(targetLine, exactFlux);
    }

    {
        std::vector<double> targetLine{1., 0., 0.,
                                       EPS, 1., 0.};
        double exactFlux = 11;
        test1Triangle(targetLine, exactFlux);
    }
    
    {
        std::vector<double> targetLine{0., 1., 0.,
                                       0., 0., 0.};
        double exactFlux = -13;
        test1Triangle(targetLine, exactFlux);
    }
    
    return 0;
}
