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

    std::cerr << "building...\n";
    int num_cells_per_bucket = 1;
    ier = mnt_regridedges2d_build(&rg, num_cells_per_bucket);
    assert(ier == 0);
    std::cerr << "done!\n";

    ier = mnt_regridedges2d_print(&rg);
    assert(ier == 0);    
    
    ier = mnt_regridedges2d_dumpWeights(&rg, weightsFileName.c_str(), weightsFileName.size());
    assert(ier == 0);

    ier = mnt_regridedges2d_del(&rg);
    assert(ier == 0);

}

int main() {

    /*

    tiny1x1 -> tiny1x1

          2
    3 ----<-----2
    |           |
 3  ^           V 1
    |           |
    0 ---->-----1
          0
    */

    //test("@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics", 
    //     "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics",
    //     "weightsTiny1x1_to_tiny1x1.nc");

    /*
    tiny1x1 -> tiny1x2


          2          6
    3 ----<-----2 --->-----5
    |           |          |
 3  ^           V 1        ^  5   
    |           |          |
    0 ---->-----1---->-----4
          0          4

    */

    test("@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics", 
         "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc:physics",
         "weightsTiny1x1_to_tiny1x2.nc");


}