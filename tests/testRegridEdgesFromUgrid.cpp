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

void test2() {

    int ier;
    std::string srcFile = "@CMAKE_SOURCE_DIR@/data/cs_16.nc";
    std::string dstFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";
    std::string outputFile = "weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrc(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << "test2: loadSrc...OK\n";

    ier = mnt_regridedges_loadDst(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << "test2: loadDst...OK\n";

    int numCellsPerBucket = 8;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);
    std::cerr << "test2: build...OK\n";

    vtkIdType edgeId;
    int edgeSign;
    double p0[3];
    double p1[3];

    // set the source data
    ier = mnt_grid_new(&rg->srcGridObj);
    ier = mnt_grid_loadFrom2DUgrid(&rg->srcGridObj, srcFile.c_str());
    size_t numSrcCells, numSrcEdges;
    ier = mnt_grid_getNumberOfCells(&rg->srcGridObj, &numSrcCells);
    ier = mnt_grid_getNumberOfUniqueEdges(&rg->srcGridObj, &numSrcEdges);
    std::vector<double> srcDataExact(numSrcEdges);
    for (size_t cellId = 0; cellId < numSrcCells; ++cellId) {
        for (int ie = 0; ie < 4; ++ie) {
            ier = mnt_grid_getPoints(&rg->srcGridObj, cellId, ie, p0, p1);
            ier = mnt_grid_getEdgeId(&rg->srcGridObj, cellId, ie, &edgeId, &edgeSign);
            srcDataExact[edgeId] = edgeSign * (streamFunc(p1) - streamFunc(p0));
        }
    }

    // set the exact destination data
    ier = mnt_grid_new(&rg->dstGridObj);
    ier = mnt_grid_loadFrom2DUgrid(&rg->dstGridObj, dstFile.c_str());
    size_t numDstCells, numDstEdges;
    ier = mnt_grid_getNumberOfCells(&rg->dstGridObj, &numDstCells);
    ier = mnt_grid_getNumberOfUniqueEdges(&rg->dstGridObj, &numDstEdges);
    std::vector<double> dstDataExact(numDstEdges);
    for (size_t cellId = 0; cellId < numDstCells; ++cellId) {
        for (int ie = 0; ie < 4; ++ie) {
            ier = mnt_grid_getPoints(&rg->dstGridObj, cellId, ie, p0, p1);
            ier = mnt_grid_getEdgeId(&rg->dstGridObj, cellId, ie, &edgeId, &edgeSign);
            dstDataExact[edgeId] = edgeSign * (streamFunc(p1) - streamFunc(p0));
        }
    }

    // regrid
    std::vector<double> dstData(numDstEdges);
    ier = mnt_regridedges_applyWeightsToEdgeIdField(&rg, &srcDataExact[0], numDstEdges, &dstData[0]);
    assert(ier == 0);

    // check
    double totError = 0;
    for (size_t i = 0; i < numDstEdges; ++i) {
        std::cerr << " edge Id: " << i << " dstData = " << dstData[i] << " exact: " << dstDataExact[i] << '\n';
        totError += std::abs(dstData[i] - dstDataExact[i]);
    }

    std::cout << " total interpolation |error|: " << totError << '\n';

    // clean up
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

}


int main() {

    test1();
    test2();

    return 0;
}   
