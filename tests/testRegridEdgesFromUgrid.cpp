#include <mntRegridEdges.h>
#undef NDEBUG // turn on asserts

double streamFunc(const double p[]) {
    return p[0];
}

void test1() {

    int ier;
    std::string srcFile = "@CMAKE_SOURCE_DIR@/data/cs_16.nc";
    std::string dstFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";
    std::string outputFile = "weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrc(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_loadDst(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);

    int numCellsPerBucket = 8;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);

    ier = mnt_regridedges_dump(&rg, outputFile.c_str(), (int) outputFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

}

void regridTest(const std::string& testName, const std::string& srcFile, const std::string& dstFile) {

    int ier;
    std::string outputFile = testName + "Weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrc(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadSrc...OK\n";

    ier = mnt_regridedges_loadDst(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadDst...OK\n";

    int numCellsPerBucket = 1;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);
    std::cerr << testName << ": build...OK\n";

    vtkIdType edgeId;
    int edgeSign;
    double p0[3];
    double p1[3];

    // set the source data
    size_t numSrcCells;
    ier = mnt_grid_getNumberOfCells(&rg->srcGridObj, &numSrcCells);
    std::vector<double> srcData(numSrcCells * 4);
    for (size_t cellId = 0; cellId < numSrcCells; ++cellId) {
        for (int ie = 0; ie < 4; ++ie) {
            ier = mnt_grid_getPoints(&rg->srcGridObj, cellId, ie, p0, p1);
            srcData[cellId*4 + ie] = streamFunc(p1) - streamFunc(p0);
        }
    }

    size_t numDstCells;
    ier = mnt_grid_getNumberOfCells(&rg->dstGridObj, &numDstCells);
    std::vector<double> dstData(numDstCells * 4);

    // regrid
    ier = mnt_regridedges_applyWeightsToCellEdgeField(&rg, &srcData[0], &dstData[0]);
    assert(ier == 0);

    // check
    double totError = 0;
    for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
        for (int ie = 0; ie < 4; ++ie) {
            size_t k = dstCellId*4 + ie;
            std::cerr << testName << ": dst edge Id: " << dstCellId << " edge " << ie << " dstData = " << dstData[k] << " exact: " << srcData[k] << '\n';
            totError += std::abs(dstData[k] - srcData[k]);
        }
    }

    std::cout << testName << ": total interpolation |error|: " << totError << '\n';
    assert(totError < 1.e-12);

    // clean up
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);


}


int main() {

    test1();
    regridTest("tiny1x2_1x1", "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc", "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc");
    regridTest("tiny1x1_1x2", "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc", "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc");
    regridTest("same", "@CMAKE_SOURCE_DIR@/data/cs_4.nc", "@CMAKE_SOURCE_DIR@/data/cs_4.nc"); 
    regridTest("cs16_4", "@CMAKE_SOURCE_DIR@/data/cs_16.nc", "@CMAKE_SOURCE_DIR@/data/cs_4.nc"); 

    return 0;
}   
