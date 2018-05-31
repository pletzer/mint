#include <mntRegridEdges.h>
#include <iostream>
#include <vtkIdList.h>
#include <MvVector.h>

extern "C"
int mnt_regridedges_new(mntLatLon_t** self) {
    *self = new mntRegridEdges_t();
    (*self)->srcGrid = vtkUnstructuredGrid::New();
    (*self)->dstGrid = vtkUnstructuredGrid::New();
    (*self)->srcLoc = vtkCellLocator::New();
    return 0;
}

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(mntLatLon_t** self) {
    (*self)->srcGrid = vtkUnstructuredGrid::Delete();
    (*self)->dstGrid = vtkUnstructuredGrid::Delete();
    (*self)->srcLoc = vtkCellLocator::Delete();
    *self = new mntRegridEdges_t();
    return 0;
}

extern "C"
int mnt_regridedges_setSrcGrid(mntLatLon_t** self, vtkUnstructuredGrid* grid) {
    (*self)->srcGrid = srcGrid;
}

extern "C"
int mnt_regridedges_setDstGrid(mntLatLon_t** self, vtkUnstructuredGrid* grid) {
    (*self)->srcGrid = srcGrid;
}

extern "C"
int mnt_regridedges_build(mntLatLon_t** self) {
    (*self)->srcLoc.SetDataSet((*self)->srcGrid);
    (*self)->srcLoc.BuildLocator(); 

    // compute the weights
    vtkIdList* dstPtIds = vtkIdList::New();
    vtkIdList* srcCellIds = vtkIdList::New();
    double* dstEdgePt0;
    double* dstEdgePt1;

    vtkIdType numSrcCells = (*self)->srcGrid.GetNumberOfCells();
    vtkIdType numDstCells = (*self)->dstGrid.GetNumberOfCells();

    // iterate over the dst grid cells
    for (vtkIdType dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {

        // iterate over the four edges of each dst cell
        (*self)->dstGrid.GetCellPoints(dstCellId, dstPtIds);
        for (size_t i0 = 0; i0 < 4; ++i0) {

            size_t i1 = (i0 + 1) % 4;

            // get the start/end points of the dst edge
            vtkIdType id0 = dstPtIds.GetId(i0);
            vtkIdType id1 = dstPtIds.GetId(i1);
                
            dstEdgePt0 = (*self)->dstGrid.GetPoint(id0);
            dstEdgePt1 = (*self)->dstGrid.GetPoint(id1);

            // break the edge into sub-edges
            mntPolysegmentIter bsi = mntPolysegmentIter((*self)->srcGrid, 
                                                        (*self)->srcLoc,
                                                        dstEdgePt0, dstEdgePt1);

            size_t numSegs = bsi.getNumSegs();
            for (size_t iseg = 0; iseg < numSegs; ++iseg) {

                vtkIdType srcCellId = bsi.getCellId();
                const MvVector<double>& xia = bsi.getBegCellParamCoord();
                const MvVector<double>& xib = seg.getEndCellParamCoord();
                const double& coeff = bsi.getCoefficient();

                MvVector<double> dxi = xib - xia;
                MvVector<double> xiMid = 0.5*(xia + xib);

                std::pair<vtkIdType, vtkIdType> k = std::pair(dstCellId, srcCellId);

                // compute the weight contributions from each src edge
                double ws[] = {+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                               + dxi[1] * (0.0 + xiMid[0]) * coeff,
                               - dxi[0] * (0.0 + xiMid[1]) * coeff,
                               - dxi[1] * (1.0 - xiMid[0]) * coeff};

                if ((*self)->weights.find(k) == (*self)->weights.end()) {
                    // initialize the weights to z zero 4x4 matrix
                    ColMat<double> zeros(4, 4, 0.0);
                    std::pair< std::pair<vtkIdType, vtkIdType> > kv = std::pair(k, zeros);
                    (*self)->weights.insert(kv);
                }

                std::map< std::pair<vtkIdType, vtkIdType>, ColMat<double> >::iterator 
                     it = (*self)->weights.find(k);

                // add the weights
                for (size_t j = 0; j < it.second.size(1); ++j)
                    it.second(i0, j) += ws[j];

                bsi.next();

            }

            double totalT = bsi.getIntegratedParamCoord();
            if abs(totalT - 1.0) > 1.e-6:
                    print('Warning: total t of segment: {:.3f} != 1 (diff={:.1g}), dst cell {} points=[{}, {}]'.format(\
                        totalT, totalT - 1.0, dstCellId, dstEdgePt0, dstEdgePt1))

        }
    }

    // clean up
    srcCellIds->Delete();
    dstPtIds->Delete();

}

/**
 * Load the weights from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_load(mntLatLon_t** self, const std::string& filename);

/**
 * Dump the weights to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dump(mntLatLon_t** self, const std::string& filename);


