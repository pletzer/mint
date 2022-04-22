#include "mntRegridEdges.h"
#include "mntLogger.h"
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>

#define EPS 0.01

void test2Triangles() {
    
    int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    /* 2 triangles

     5<2.4.....<6......3
     :    \            :
     v     \           ^
     3      \          5
     :       ^         :
     :        1        :
     :         \       :
     :          \      :
     :           \     :
     :            \    :
     0.......0>....1.4>2
    */

    std::size_t ncells = 2;
    std::size_t nedges = 7;
    std::size_t npoints = 6;

    std::vector<double> xyz{0.    ,0.,0.,
                            1.-EPS,0.,0.,
                            1.    ,0.,0.,
                            1.    ,1.,0.,
                            EPS   ,1.,0.,
                            0.    ,1.,0.};

    std::vector<std::size_t> face2nodes{0, 1, 4, 5,
                                        2, 3, 4, 1};

    std::vector<std::size_t> edge2nodes{0, 1,
                                        1, 4,
                                        4, 5,
                                        5, 0,
                                        1, 2,
                                        2, 3,
                                        3, 4};
    ier = mnt_grid_loadFromUgrid2DData(&grd, ncells, nedges, npoints, &xyz[0], &face2nodes[0], &edge2nodes[0]);
    assert(ier == 0);

    RegridEdges_t* rg = NULL;
    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    // src and dst grids are the same
    ier = mnt_regridedges_setSrcGrid(&rg, grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rg, grd);
    assert(ier == 0);

    int num_cells_per_bucket = 128;
    double periodX = 0.;
    int enableFolding = 0;
    ier = mnt_regridedges_buildLocator(&rg, num_cells_per_bucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rg, debug);
    assert(ier == 0);

    // data values are the edge indices
    std::size_t numEdges;
    ier = mnt_regridedges_getNumSrcEdges(&rg, &numEdges);
    assert(ier == 0);

    // set the edge field integrals to the edge Ids
    std::vector<double> srcData(numEdges);
    std::vector<double> dstData(numEdges);
    for (auto i = 0; i < numEdges; ++i) {
        srcData[i] = i;
    }
    // set the fluxes to zero on the small edges
    srcData[2] = 0;
    srcData[4] = 0;

    // apply the weights
    int placement = MNT_UNIQUE_EDGE_DATA;
    ier = mnt_regridedges_apply(&rg, &srcData[0], &dstData[0], placement);
    assert(ier == 0);

    mnt_printLogMessages();

    // check
    double error = 0;
    for (auto i = 0; i < numEdges; ++i) {
        error += std::abs(dstData[i] - srcData[i]);
        std::cout << i << " srcData = " << srcData[i] << " dstData = " << dstData[i] << '\n';
    }
    error /= (double) numEdges;

    std::cerr << "2 triangles error: " << error << '\n';
    assert(error < 1.e-10);
   
    // reclaim the memory
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);
}

void test1Quad() {
    
    int ier;
    Grid_t* grd;
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
    assert(ier == 0);

    RegridEdges_t* rg = NULL;
    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    // lat-lon
    int fixLonAcrossDateline = 0;
    int averageLonAtPole = 0;
    ier = mnt_regridedges_setSrcGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);
    ier = mnt_regridedges_setDstGridFlags(&rg, fixLonAcrossDateline, averageLonAtPole);
    assert(ier == 0);

    // src and dst grids are the same
    ier = mnt_regridedges_setSrcGrid(&rg, grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rg, grd);
    assert(ier == 0);

    int num_cells_per_bucket = 128;
    double periodX = 360.;
    int enableFolding = 0;
    ier = mnt_regridedges_buildLocator(&rg, num_cells_per_bucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rg, debug);
    assert(ier == 0);

    // data values are the edge indices
    std::size_t numEdges;
    ier = mnt_regridedges_getNumSrcEdges(&rg, &numEdges);
    assert(ier == 0);

    std::vector<double> srcData(numEdges);
    std::vector<double> dstData(numEdges);
    for (auto i = 0; i < numEdges; ++i) {
        srcData[i] = i;
    }

    // apply the weights
    int placement = MNT_UNIQUE_EDGE_DATA;
    ier = mnt_regridedges_apply(&rg, &srcData[0], &dstData[0], placement);
    assert(ier == 0);

    mnt_printLogMessages();

    // check
    double error = 0;
    for (auto i = 0; i < numEdges; ++i) {
        error += std::abs(srcData[i] - dstData[i]);
    }
    error /= (double) numEdges;

    std::cerr << "lat-lon error: " << error << '\n';
    assert(error < 1.e-10);
   
    // reclaim the memory
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);
}

int main() {
    test2Triangles();
    test1Quad();
}
