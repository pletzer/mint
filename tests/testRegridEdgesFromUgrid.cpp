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

void testSameGrid4() {

    int ier;
    std::string srcFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";
    std::string dstFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";
    std::string outputFile = "weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrc(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << "testSameGrid4: loadSrc...OK\n";

    ier = mnt_regridedges_loadDst(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << "testSameGrid4: loadDst...OK\n";

    int numCellsPerBucket = 1;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);
    std::cerr << "testSameGrid4: build...OK\n";

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
            std::cerr << " dst edge Id: " << dstCellId << " edge " << ie << " dstData = " << dstData[k] << " exact: " << srcData[k] << '\n';
            totError += std::abs(dstData[k] - srcData[k]);
        }
    }

    std::cout << " total interpolation |error|: " << totError << '\n';
    assert(totError < 1.e-12);

    // clean up
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
    testSameGrid4();
    test2();

    return 0;
}   
