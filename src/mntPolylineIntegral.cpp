#include "mntLogger.h"
#include <mntGlobal.h>
#include <mntPolylineIntegral.h>
#include <mntPolysegmentIter.h>
#include <mntWeights.h>
#include <iostream>

//#define DEBUG

LIBRARY_API
int mnt_polylineintegral_new(PolylineIntegral_t** self) {
    
    *self = new PolylineIntegral_t();

    (*self)->loc = vmtCellLocator::New();
    (*self)->vgrid = NULL;

    return 0;
}

LIBRARY_API
int mnt_polylineintegral_del(PolylineIntegral_t** self) {

    // destroy everything here
    (*self)->loc->Delete();

    delete *self;

    return 0;
}


LIBRARY_API
int mnt_polylineintegral_setGrid(PolylineIntegral_t** self, Grid_t* grid) {

    int ier = 0;
    std::string msg;

    // get VTK grid
    mnt_grid_get(&grid, &(*self)->vgrid);
    (*self)->gridObj = grid;
    (*self)->loc->SetDataSet((*self)->vgrid);

    return ier;
}


LIBRARY_API
int mnt_polylineintegral_buildLocator(PolylineIntegral_t** self,
                                      int numCellsPerBucket,
                                      double periodX,
                                      int enableFolding) {

    int ier = 0;
    std::string msg;

    // build the cell locator
    (*self)->loc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    (*self)->loc->setPeriodicityLengthX(periodX);
    if (enableFolding == 1) {
        (*self)->loc->enableFolding();
    }
    (*self)->loc->BuildLocator();

    return ier;
}

LIBRARY_API
int mnt_polylineintegral_computeWeights(PolylineIntegral_t** self,
                                        int npoints, const double xyz[], int counterclock) {

    int ier = 0;
    std::string msg;

    // reset
    (*self)->weights.clear();

    if (npoints < 2) {
        std::string mgs = "need at least two points";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    if (!(*self)->vgrid) {
        std::string mgs = "need to build locator before computing the weights";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    const double* xi0_xi1 = QUAD_POSITIVEXI_EDGE_DIRECTION;
    if (counterclock > 0) {
        xi0_xi1 = QUAD_COUNTERCLOCKWISE_EDGE_DIRECTION;
    }

    // iterate over the Polyline segments
    for (int ip0 = 0; ip0 < npoints - 1; ++ip0) {

        // get the starting point
        std::size_t k0 = 3*ip0;
        double p0[] = {xyz[k0 + 0], xyz[k0 + 1], xyz[k0 + 2]};

        // get the end point
        std::size_t k1 = k0 + 3;        
        double p1[] = {xyz[k1 + 0], xyz[k1 + 1], xyz[k1 + 2]};

        // build the polysegment iterator. This breaks segment p0 -> p1 into 
        // smaller segments, each being entirely contained within grid cells
        PolysegmentIter polyseg((*self)->vgrid, (*self)->loc, p0, p1);

        // number of sub-segments
        std::size_t numSubSegs = polyseg.getNumberOfSegments();
        polyseg.reset();

        // iterate over the sub-segments 
        for (std::size_t i = 0; i < numSubSegs; ++i) {

            // get the cell Id to which this sub-segment belongs to
            vtkIdType cellId = polyseg.getCellId();

            // get the cell parametric start/end values for the sub-segment
            const Vec3& xia = polyseg.getBegCellParamCoord();
            const Vec3& xib = polyseg.getEndCellParamCoord();

            // some sub-segments belong to two or more cells (e.g. if the sub-segment
            // runs along a cell edge)
            double coeff = polyseg.getCoefficient();

            // store the weights for this sub-segment
            for (int edgeIndex = 0; edgeIndex < MNT_NUM_EDGES_PER_QUAD; ++edgeIndex) {

                // starting point of the cell edge in parameter space
                const double* xi0 = &xi0_xi1[6*edgeIndex];

                // end point of the cell edge in parameter space
                const double* xi1 = &xi0_xi1[6*edgeIndex + 3];

                double weight = computeWeight(xi0, xi1, xia, xib);
#ifdef DEBUG
                msg = "ip0: " + std::to_string(ip0) + " sgmt: " + std::to_string(i) + " cell " + std::to_string(cellId) + " edge: " + std::to_string(edgeIndex) + 
                      " xi beg/end: " + std::to_string(xia[0]) + "," + std::to_string(xia[1]) + ' ' + 
                      std::to_string(xib[0]) + "," + std::to_string(xib[1]) + " wght: " + std::to_string(weight);
                mntlog::info(__FILE__, __func__, __LINE__, msg);
#endif
                // increment the weight. note the default value is 0 for numeric values
                std::pair<vtkIdType, int> ce(cellId, edgeIndex);
                (*self)->weights[ce] += weight * coeff;
            }

            // increment the iterator
            polyseg.next();       
        }

        // check if the line segment is fully accounted for
        double tTotal = polyseg.getIntegratedParamCoord();
        const double tol = 1.e-10;
        if (std::abs(tTotal - 1.0) > tol) {
            msg = "total integrated length for segment " + std::to_string(ip0) + " is " + 
                  std::to_string(tTotal) + " != 1 (diff=" + std::to_string(tTotal - 1.) + ")";
            mntlog::warn(__FILE__, __func__, __LINE__, msg);
            ier++;
        }

    }

    return ier;
}

LIBRARY_API
int mnt_polylineintegral_getIntegralOfCellByCellData(PolylineIntegral_t** self, const double data[],
                                                     double* result) {

    *result = 0;
    for (auto it = (*self)->weights.cbegin(); it != (*self)->weights.cend(); ++it) {

        // get the cell Id
        vtkIdType cellId = it->first.first;

        // get the edge index
        int edgeIndex = it->first.second;

        // get the weight
        double wght = it->second;

        // add the contribution
#ifdef DEBUG
        std::string msg = "adding weight = " + std::to_string(wght) + " * edge integral = " + 
                          std::to_string(data[cellId*MNT_NUM_EDGES_PER_QUAD + edgeIndex]) + " to the flux (" +
                          std::to_string(*result) + " so far) for cellId = " + 
                          std::to_string(cellId) + " edgeIndex = " + std::to_string(edgeIndex);
        mntlog::info(__FILE__, __func__, __LINE__, msg);
#endif
        *result += wght * data[cellId*MNT_NUM_EDGES_PER_QUAD + edgeIndex];
    }
#ifdef DEBUG
    std::string msg = "total flux = " + std::to_string(*result);
    mntlog::info(__FILE__, __func__, __LINE__, msg);
#endif

    return 0;
}

LIBRARY_API
int mnt_polylineintegral_getIntegralOfUniqueEdgeData(PolylineIntegral_t** self, const double data[],
                                                     double* result) {

    *result = 0;
    int ier = 0;

    // make sure (*self)->srcGridObj.faceNodeConnectivity and the rest have been allocated
    if (!(*self)->gridObj ||
        (*self)->gridObj->faceNodeConnectivity.size() == 0 || 
        (*self)->gridObj->faceEdgeConnectivity.size() == 0 ||
        (*self)->gridObj->edgeNodeConnectivity.size() == 0) {
        std::string msg = "grid connectivity not set (?) May need to read the grid from a netcdf Ufile";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    for (auto it = (*self)->weights.cbegin(); it != (*self)->weights.cend(); ++it) {

        // get the cell Id
        vtkIdType cellId = it->first.first;
        // get the edge index
        int edgeIndex = it->first.second;

        std::size_t edgeId;
        int edgeSign;
        ier = mnt_grid_getEdgeId(&((*self)->gridObj), cellId, edgeIndex, &edgeId, &edgeSign);

        // get the weight
        double wght = it->second;

        // add the contribution
#ifdef DEBUG
        std::string msg = "adding weight = " + std::to_string(wght) + " * edge integral = " + 
                          std::to_string(data[edgeId] * edgeSign) + " to the flux (" +
                          std::to_string(*result) + " so far) for cellId = " + 
                          std::to_string(cellId) + " edgeIndex = " + std::to_string(edgeIndex);
        mntlog::info(__FILE__, __func__, __LINE__, msg);
#endif
        *result += wght * data[edgeId] * edgeSign;
    }
#ifdef DEBUG
    std::string msg = "total flux = " + std::to_string(*result);
    mntlog::info(__FILE__, __func__, __LINE__, msg);
#endif

    return ier;
}

LIBRARY_API
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self, const double data[],
                                     int placement, double* result) {

    if (placement == MNT_CELL_BY_CELL_DATA) {
        return mnt_polylineintegral_getIntegralOfCellByCellData(self, data, result);
    }
    else {
        return mnt_polylineintegral_getIntegralOfUniqueEdgeData(self, data, result);
    }
}
