#include <mntRegridEdges.h>
#undef NDEBUG // turn on asserts

int main() {
    int ier;
    std::string srcFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";
    std::string dstFile = "@CMAKE_SOURCE_DIR@/data/cs_16.nc";
    std::string outputFile = "weights.nc";
    std::cerr << "0\n";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);
    std::cerr << "1\n";

    ier = mnt_regridedges_loadSrc(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << "2\n";

    ier = mnt_regridedges_loadDst(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << "3\n";

    int numCellsPerBucket = 8;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);
    std::cerr << "4\n";

    ier = mnt_regridedges_dump(&rg, outputFile.c_str(), (int) outputFile.size());
    assert(ier == 0);
    std::cerr << "5\n";

    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);
    std::cerr << "6\n";

    return 0;
}   
