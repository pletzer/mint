#include "mntRegridEdges.h"
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>

void testIdentity() {
    
    std::string fileMesh = "@CMAKE_SOURCE_DIR@/data/latlon4x2.nc$latlon";

    RegridEdges_t* rg = NULL;
    int ier;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    // src and dst grids are the same

    ier = mnt_regridedges_loadSrcGrid(&rg, fileMesh.c_str(), (int) fileMesh.size());
    assert(ier == 0);
    ier = mnt_regridedges_loadDstGrid(&rg, fileMesh.c_str(), (int) fileMesh.size());
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

    // check
    double error = 0;
    for (auto i = 0; i < numEdges; ++i) {
        error += std::abs(srcData[i] - dstData[i]);
    }
    std::cerr << "Error: " << error << '\n';
    assert(error < 1.e-10);
   
    // reclaim the memory
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);
}


int main() {
    testIdentity();
}
