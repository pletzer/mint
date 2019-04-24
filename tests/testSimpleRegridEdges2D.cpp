#include "mntRegridEdges2D.h"
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>


void test(const std::string& srcGridFileMesh, 
          const std::string& dstGridFileMesh,
          const std::string& weightsFileName) {

    RegridEdges2D_t* rg;
    int ier;

    ier = mnt_regridedges2d_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges2d_loadSrcGrid(&rg, srcGridFileMesh.c_str(), srcGridFileMesh.size());
    assert(ier == 0);

    ier = mnt_regridedges2d_loadDstGrid(&rg, dstGridFileMesh.c_str(), dstGridFileMesh.size());
    assert(ier == 0);

    int num_cells_per_bucket = 1;
    ier = mnt_regridedges2d_build(&rg, num_cells_per_bucket);
    assert(ier == 0);

    ier = mnt_regridedges2d_print(&rg);
    assert(ier == 0);    
    
    ier = mnt_regridedges2d_dumpWeights(&rg, weightsFileName.c_str(), weightsFileName.size());
    assert(ier == 0);

    ier = mnt_regridedges2d_del(&rg);
    assert(ier == 0);

}

int main() {

    test("@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics", 
         "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics",
         "weightsTiny1x1_1x1.nc");

}