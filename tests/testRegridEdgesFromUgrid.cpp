#include <mntRegridEdges.h>
#include <cmath>
#undef NDEBUG // turn on asserts

double streamFunc(const double p[]) {
    double lon = p[0];
    double lat = p[1];
    return cos(M_PI*lat/180.0) * sin(M_PI*lon/180.);
}

void test1() {

    int ier;
    std::string srcFile = "@CMAKE_SOURCE_DIR@/data/vel_cs_16.nc:physics";
    std::string dstFile = "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics";
    std::string outputFile = "weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);

    size_t numSrcEdges, numDstEdges;

    ier = mnt_regridedges_getNumSrcUniqueEdges(&rg, &numSrcEdges);
    assert(ier == 0);

    ier = mnt_regridedges_getNumDstUniqueEdges(&rg, &numDstEdges);
    assert(ier == 0);

    int numCellsPerBucket = 8;
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    assert(ier == 0);

    ier = mnt_regridedges_dumpWeights(&rg, outputFile.c_str(), (int) outputFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

    // read the weights and interpolate

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);

    ier = mnt_regridedges_loadWeights(&rg, outputFile.c_str(), (int) outputFile.size());
    assert(ier == 0);

    std::vector<double> srcData(numSrcEdges);
    std::vector<double> dstData(numDstEdges);

    // load the source data from file
    std::string fieldName = "line_integrated_velocity";
    ier = mnt_regridedges_loadUniqueEdgeField(&rg, srcFile.c_str(), (int) srcFile.size(),
                                              fieldName.c_str(), (int) fieldName.size(),
                                              numSrcEdges, &srcData[0]);
    assert(ier == 0);

    ier = mnt_regridedges_applyUniqueEdge(&rg, &srcData[0], &dstData[0]);
    assert(ier == 0);

    std::string resFile = "regridded_line_integrated_velocity.nc";
    ier = mnt_regridedges_dumpUniqueEdgeField(&rg, resFile.c_str(), (int) resFile.size(),
                                              fieldName.c_str(), (int) fieldName.size(),
                                              numDstEdges, &dstData[0]);
    assert(ier == 0);

    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

}

void regridCellEdgeFieldTest(const std::string& testName, const std::string& srcFile, const std::string& dstFile) {

    int ier;
    std::string outputFile = testName + "Weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadSrc...OK\n";

    ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadDst...OK\n";

    int numCellsPerBucket = 8;
    
    assert(ier == 0);
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    std::cerr << testName << ": build...OK\n";

    std::string weightFile = testName + "Weights.nc";
    ier = mnt_regridedges_dumpWeights(&rg, weightFile.c_str(), (int) weightFile.size());
    assert(ier == 0);


    vtkIdType edgeId;
    int edgeSign;
    double p0[3];
    double p1[3];

    // set the source data
    size_t numSrcCells;
    ier = mnt_grid_getNumberOfCells(&rg->srcGridObj, &numSrcCells);
    std::vector<double> srcData(numSrcCells * 4);
    for (size_t srcCellId = 0; srcCellId < numSrcCells; ++srcCellId) {
        for (int ie = 0; ie < 4; ++ie) {
            ier = mnt_grid_getPoints(&rg->srcGridObj, srcCellId, ie, p0, p1);

            size_t k = srcCellId*4 + ie;
            srcData[k] = streamFunc(p1) - streamFunc(p0);
        }
    }

    size_t numDstCells;
    ier = mnt_grid_getNumberOfCells(&rg->dstGridObj, &numDstCells);
    std::vector<double> dstData(numDstCells * 4);

    // regrid
    ier = mnt_regridedges_applyCellEdge(&rg, &srcData[0], &dstData[0]);
    assert(ier == 0);

    double totalLoopIntegrals = 0.0;

    // check
    printf("%s\n dstCellId         edge        interpVal      exact        error               p0               p1\n", testName.c_str());
    double totError = 0;
    for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
        
        double cellLoopIntegral = 0.0;

        for (int ie = 0; ie < 4; ++ie) {

            // get the coordinates of the start (p0) and end (p1) points. This may involve
            // adding/subtracting 360 deg to the longitude of the cell crosses the dateline
         	int circSign;
        	ier = mnt_regridedges_getDstEdgePoints(&rg, dstCellId, ie, &circSign, p0, p1);

            // exact line integral
            double exact = streamFunc(p1) - streamFunc(p0);

            // interpolated line integral 
            size_t k = dstCellId*4 + ie;
            double interpVal = dstData[k];

            // orientation for loop integral is counterclockwise
            // the first two edges are the direction of the contour
            // integral, the last two are in opposite direction
            // the value circSign relfects this.
            cellLoopIntegral += interpVal * circSign;

            // numerical error
            double error = interpVal - exact;
            if (std::abs(error) > 1.e-8) {
                printf("%10ld           %1d        %10.6lf   %10.6lf    %12.5lg     %5.1lf,%5.1lf      %5.1lf,%5.1lf\n", 
                    dstCellId, ie, interpVal, exact, error, p0[0], p0[1], p1[0], p1[1]);
            }
            totError += std::abs(error);
        }

        if (std::abs(cellLoopIntegral) > 1.e-8) {
          printf("%10ld loop integral = %12.5lg\n", dstCellId, cellLoopIntegral);
        }

        totalLoopIntegrals += std::abs(cellLoopIntegral);
    }

    std::cout << testName << ": total interpolation |error|: " << totError << '\n';
    std::cout << testName << ": sum of |loop integrals|    : " << totalLoopIntegrals << '\n';
    assert(totError < 1.e-8);
    assert(totalLoopIntegrals  < 1.e-8);

    // clean up
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);


}

void regridUniqueEdgeIdFieldTest(const std::string& testName, const std::string& srcFile, const std::string& dstFile) {

    int ier;
    std::string outputFile = testName + "Weights.nc";

    RegridEdges_t* rg;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_loadSrcGrid(&rg, srcFile.c_str(), (int) srcFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadSrc...OK\n";

    ier = mnt_regridedges_loadDstGrid(&rg, dstFile.c_str(), (int) dstFile.size());
    assert(ier == 0);
    std::cerr << testName << ": loadDst...OK\n";

    int numCellsPerBucket = 8;
    
    assert(ier == 0);
    ier = mnt_regridedges_build(&rg, numCellsPerBucket);
    std::cerr << testName << ": build...OK\n";

    std::string weightFile = testName + "Weights.nc";
    ier = mnt_regridedges_dumpWeights(&rg, weightFile.c_str(), (int) weightFile.size());
    assert(ier == 0);


    vtkIdType edgeId;
    int edgeSign;
    double p0[3];
    double p1[3];

    // set the source data
    size_t numSrcCells, numSrcEdges;
    vtkIdType srcEdgeId;
    int srcEdgeSign;
    ier = mnt_grid_getNumberOfCells(&rg->srcGridObj, &numSrcCells);
    assert(ier == 0);
    ier = mnt_grid_getNumberOfUniqueEdges(&rg->srcGridObj, &numSrcEdges);
    assert(ier == 0);

    std::vector<double> srcData(numSrcEdges);

    for (size_t srcCellId = 0; srcCellId < numSrcCells; ++srcCellId) {
        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getEdgeId(&rg->srcGridObj, srcCellId, ie, &srcEdgeId, &srcEdgeSign);
            assert(ier == 0);

            ier = mnt_grid_getPoints(&rg->srcGridObj, srcCellId, ie, p0, p1);
            assert(ier == 0);

            srcData[srcEdgeId] = srcEdgeSign * (streamFunc(p1) - streamFunc(p0));
        }
    }

    size_t numDstCells, numDstEdges;
    vtkIdType dstEdgeId;
    int dstEdgeSign;
    ier = mnt_grid_getNumberOfCells(&rg->dstGridObj, &numDstCells);
    assert(ier == 0);
    ier = mnt_grid_getNumberOfUniqueEdges(&rg->dstGridObj, &numDstEdges);
    assert(ier == 0);

    std::vector<double> dstData(numDstEdges);

    // regrid
    ier = mnt_regridedges_applyUniqueEdge(&rg, &srcData[0], &dstData[0]);
    assert(ier == 0);

    // check
    printf("%s\n dstCellId         edgeIndex        edgeId      interpVal      exact        error               p0               p1\n", testName.c_str());
    double totError = 0;
    for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {
        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getEdgeId(&rg->dstGridObj, dstCellId, ie, &dstEdgeId, &dstEdgeSign);
            assert(ier == 0);

            ier = mnt_grid_getPoints(&rg->dstGridObj, dstCellId, ie, p0, p1);
            assert(ier == 0);

            double exact = dstEdgeSign * (streamFunc(p1) - streamFunc(p0));

            double interpVal = dstData[dstEdgeId];

            double error = interpVal - exact;

            if (std::abs(error) > 1.e-6) {
                printf("%10ld           %1d         %9lld      %10.6lf   %10.6lf    %12.5lg     %5.1lf,%5.1lf      %5.1lf,%5.1lf\n", 
                    dstCellId, ie, dstEdgeId, interpVal, exact, error, p0[0], p0[1], p1[0], p1[1]);
            }
            totError += std::abs(error);
        }
    }

    std::cout << testName << ": total interpolation |error|: " << totError << '\n';
    assert(totError < 1.e-8);

    // clean up
    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

}


int main() {

    test1();

    // crashes when building the cell locator
    //regridTest("tiny1x2_1x1", "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc:physics", "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics");
    //regridTest("tiny1x1_1x2", "@CMAKE_SOURCE_DIR@/data/tiny1x1.nc:physics", "@CMAKE_SOURCE_DIR@/data/tiny1x2.nc:physics");

    regridCellEdgeFieldTest("sameCellField", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics"); 
    regridUniqueEdgeIdFieldTest("sameUniqueEdgeIdField", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics");
    regridCellEdgeFieldTest("uniqueEdgeIdField_16->4", "@CMAKE_SOURCE_DIR@/data/cs_16.nc:physics", "@CMAKE_SOURCE_DIR@/data/cs_4.nc:physics"); 

    return 0;
}   
