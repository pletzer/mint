#include "mntRegridEdges2D.h"
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>



int main() {

    std::cerr << "*** 0\n";
    RegridEdges2D_t* rg;
    int ier;

    std::cerr << "*** 1\n";
    ier = mnt_regridedges2d_new(&rg);
    assert(ier == 0);

    std::cerr << "*** 2\n";
    std::string srcGridFile = "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics";
    ier = mnt_regridedges2d_loadSrcGrid(&rg, srcGridFile.c_str(), srcGridFile.size());
    assert(ier == 0);

    std::cerr << "*** 3\n";
    std::string dstGridFile = "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc:physics";
    ier = mnt_regridedges2d_loadDstGrid(&rg, dstGridFile.c_str(), dstGridFile.size());
    assert(ier == 0);

    std::cerr << "*** 4\n";
    int num_cells_per_bucket = 1;
    ier = mnt_regridedges2d_build(&rg, num_cells_per_bucket);
    assert(ier == 0);

    std::cerr << "*** 5\n";
    ier = mnt_regridedges2d_print(&rg);
    assert(ier == 0);    
    
    std::cerr << "*** 6\n";
    std::string outputFilename = "simpleRegridEdges2DWeights.nc";
    ier = mnt_regridedges2d_dumpWeights(&rg, outputFilename.c_str(), outputFilename.size());
    assert(ier == 0);

    std::cerr << "*** 7\n";
    ier = mnt_regridedges2d_del(&rg);
    assert(ier == 0);
}