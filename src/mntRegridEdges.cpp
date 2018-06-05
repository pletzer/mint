#include <mntRegridEdges.h>
#include <mntPolysegmentIter.h>
#include <iostream>
#include <vtkIdList.h>
#include <netcdf.h>

extern "C"
int mnt_regridedges_new(RegridEdges_t** self) {
    *self = new RegridEdges_t();
    (*self)->srcGrid = NULL;
    (*self)->dstGrid = NULL;
    (*self)->srcLoc = vtkCellLocator::New();
    return 0;
}

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(RegridEdges_t** self) {
    // the src and dst grids don't belong to this class
    (*self)->srcLoc->Delete();
    delete *self;
    return 0;
}

extern "C"
int mnt_regridedges_setSrcGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid) {
    // check
    if (grid->GetNumberOfPoints() != 4 * grid->GetNumberOfCells()) {
        std::cerr << "mnt_regridedges_setSrcGrid: ERROR num points != 4 * num cells!\n";
        return 1;
    }

    // borrow pointer
    (*self)->srcGrid = grid;

    return 0;
}

extern "C"
int mnt_regridedges_setDstGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid) {
    // check
    if (grid->GetNumberOfPoints() != 4 * grid->GetNumberOfCells()) {
        std::cerr << "mnt_regridedges_setDstGrid: ERROR num points != 4 * num cells!\n";
        return 1;
    }

    // borrow pointer
    (*self)->dstGrid = grid;
    return 0;
}

extern "C"
int mnt_regridedges_build(RegridEdges_t** self) {

    // checks
    if (!(*self)->srcGrid) {
        std::cerr << "mnt_regridedges_build: ERROR must set source grid!\n";
        return 1;
    }
    if (!(*self)->dstGrid) {
        std::cerr << "mnt_regridedges_build: ERROR must set destination grid!\n";
        return 2;
    }

    // build the locator
    (*self)->srcLoc->SetDataSet((*self)->srcGrid);
    (*self)->srcLoc->BuildLocator(); 

    // compute the weights
    vtkIdList* dstPtIds = vtkIdList::New();
    vtkIdList* srcCellIds = vtkIdList::New();
    double dstEdgePt0[] = {0., 0., 0.};
    double dstEdgePt1[] = {0., 0., 0.};

    vtkIdType numSrcCells = (*self)->srcGrid->GetNumberOfCells();
    vtkIdType numDstCells = (*self)->dstGrid->GetNumberOfCells();

    // iterate over the dst grid cells
    for (vtkIdType dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {

        // get this cell vertex Ids
        (*self)->dstGrid->GetCellPoints(dstCellId, dstPtIds);

        // iterate over the four edges of each dst cell
        for (size_t i0 = 0; i0 < 4; ++i0) {

            size_t i1 = (i0 + 1) % 4;

            // get the start/end points of the dst edge
            vtkIdType id0 = dstPtIds->GetId(i0);
            vtkIdType id1 = dstPtIds->GetId(i1);
                
            (*self)->dstGrid->GetPoints()->GetPoint(id0, dstEdgePt0);
            (*self)->dstGrid->GetPoints()->GetPoint(id1, dstEdgePt1);

            // break the edge into sub-edges
            PolysegmentIter polySegIter = PolysegmentIter((*self)->srcGrid, 
                                                          (*self)->srcLoc,
                                                          dstEdgePt0, dstEdgePt1);

            std::cout << "dst cell " << dstCellId
                      << " dstEdgePt0=" << dstEdgePt0[0] << ',' << dstEdgePt0[1]
                      << " dstEdgePt1=" << dstEdgePt1[0] << ',' << dstEdgePt1[1] << '\n';

            // number of sub-segments
            size_t numSegs = polySegIter.getNumberOfSegments();

            // iterate over the sub-segments. Each sub-segment gets a src cell Id,
            //start/end cell param coords, the coefficient...
            polySegIter.reset();
            for (size_t iseg = 0; iseg < numSegs; ++iseg) {

                const vtkIdType srcCellId = polySegIter.getCellId();
                const std::vector<double>& xia = polySegIter.getBegCellParamCoord();
                const std::vector<double>& xib = polySegIter.getEndCellParamCoord();
                const double coeff = polySegIter.getCoefficient();

                std::vector<double> dxi({xib[0] - xia[0], xib[1] - xia[1]});
                std::vector<double> xiMid({0.5*(xia[0] + xib[0]), 0.5*(xia[1] + xib[1])});

                std::pair<vtkIdType, vtkIdType> k = std::pair<vtkIdType, vtkIdType>(dstCellId, 
                                                                                    srcCellId);

                // compute the weight contributions from each src edge
                double ws[] = {+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                               + dxi[1] * (0.0 + xiMid[0]) * coeff,
                               - dxi[0] * (0.0 + xiMid[1]) * coeff,
                               - dxi[1] * (1.0 - xiMid[0]) * coeff};

                std::cout << "\t seg=" << iseg << " srcCellId=" << srcCellId << " xia=" << xia[0] << ',' << xia[1] << " xib=" << xib[0] << ',' << xib[1]
                          << " coeff=" << coeff << " dxi=" << dxi[0] << ',' << dxi[1] << '\n';

                if ((*self)->weights.find(k) == (*self)->weights.end()) {
                    // initialize the weights to a zero 4x4 matrix
                    std::vector<double> zeros(4, 0.0);
                    std::pair< std::pair<vtkIdType, vtkIdType>, std::vector<double> > kv 
                      = std::pair< std::pair<vtkIdType, vtkIdType>, std::vector<double> >(k, zeros);
                    (*self)->weights.insert(kv);
                }

                std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> >::iterator 
                     it = (*self)->weights.find(k);

                // add the weights
                for (size_t j = 0; j < it->second.size(); ++j) {
                    it->second[j] += ws[j];
                }

                // next segment
                polySegIter.next();

            }

            double totalT = polySegIter.getIntegratedParamCoord();
            if (std::abs(totalT - 1.0) > 1.e-6) {
                std::cerr << "Warning: total t of segment: " << totalT 
                          << " != 1 (diff=" << totalT - 1.0 
                          << " dst cell " << dstCellId << " dstEdgePt0=" 
                          << dstEdgePt0[0] << ',' << dstEdgePt0[1] << " dstEdgePt1=" 
                          << dstEdgePt1[0] << ',' << dstEdgePt1[1] << '\n';
            }

        }
    }

    // clean up
    srcCellIds->Delete();
    dstPtIds->Delete();

    return 0;
}

/**
 * Load the weights from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_load(RegridEdges_t** self, const char* filename);

/**
 * Dump the weights to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dump(RegridEdges_t** self, const char* filename) {

    size_t numWeights = (*self)->weights.size();
    const int numEdges = 4;

    int ncid, ier;
    ier = nc_create(filename, NC_CLOBBER, &ncid);

    // create dimensions
    int numWeightsId;
    int numEdgesId;
    ier = nc_def_dim(ncid, "num_weights", (int) numWeights, &numWeightsId);
    ier = nc_def_dim(ncid, "num_edges", numEdges, &numEdgesId);


    // create variables
    int numWeightsAxis[] = {numWeightsId};
    int numWeightsNumEdgesAxes[] = {numWeightsId, numEdgesId};

    int dstCellIdsId;
    ier = nc_def_var(ncid, "dst_cell_ids", NC_LONG, 1, numWeightsAxis, &dstCellIdsId);

    int srcCellIdsId;
    ier = nc_def_var(ncid, "src_cell_ids", NC_LONG, 1, numWeightsAxis, &srcCellIdsId);

    int weightsId;
    ier = nc_def_var(ncid, "weights", NC_DOUBLE, 1, numWeightsNumEdgesAxes, &weightsId);

    // close define mode
    ier = nc_enddef(ncid);

    // load into arrays
    std::vector<long> dstCellIds(numWeights);
    std::vector<long> srcCellIds(numWeights);
    std::vector<double> weights(numWeights * numEdges);
    size_t i = 0;
    for (std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> >::const_iterator 
        it = (*self)->weights.begin(); it != (*self)->weights.end(); ++it) {
        dstCellIds[i] = it->first.first;
        srcCellIds[i] = it->first.second;
        for (int j = 0; j < numEdges; ++j) {
        weights[numEdges*i + j] = it->second[j];
        }
        i++;
    }

    // write
    ier = nc_put_var_long(ncid, dstCellIdsId, &dstCellIds[0]);
    ier = nc_put_var_long(ncid, srcCellIdsId, &srcCellIds[0]);
    ier = nc_put_var_double(ncid, weightsId, &weights[0]);

    ier = nc_close(ncid);

    return 0;
}


