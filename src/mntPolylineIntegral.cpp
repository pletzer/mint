#include "mntLogger.h"
#include <mntPolylineIntegral.h>
#include <mntPolysegmentIter.h>
#include <mntWeights.h>
#include <iostream>

//#define DEBUG

LIBRARY_API
int mnt_polylineintegral_new(PolylineIntegral_t** self) {
    
    *self = new PolylineIntegral_t();

    return 0;
}

LIBRARY_API
int mnt_polylineintegral_del(PolylineIntegral_t** self) {
    int ier = 0;

    // destroy everything here

    delete *self;
    return ier;
}


LIBRARY_API
int mnt_polylineintegral_build(PolylineIntegral_t** self, Grid_t* grid, 
                               int npoints, const double xyz[], int counterclock, double periodX) {

    int ier = 0;
    std::string msg;

    if (npoints <= 0) {
        std::string mgs = "need at least one point";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    const double* xi0_xi1 = QUAD_POSITIVEXI_EDGE_DIRECTION;
    if (counterclock > 0) {
        xi0_xi1 = QUAD_COUNTERCLOCKWISE_EDGE_DIRECTION;
    }

    // get VTK grid
    vtkUnstructuredGrid* vgrid = NULL;
    mnt_grid_get(&grid, &vgrid);

    // build the cell locator
    vmtCellLocator* loc = vmtCellLocator::New();
    loc->setPeriodicityLengthX(periodX);
    loc->SetDataSet(vgrid);
    loc->BuildLocator();

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
        PolysegmentIter polyseg(vgrid, loc, p0, p1, periodX);

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
            for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex) {

                // starting point of the cell edge in parameter space
                const double* xi0 = &xi0_xi1[6*edgeIndex];

                // end point of the cell edge in parameter space
                const double* xi1 = &xi0_xi1[6*edgeIndex + 3];

                double weight = computeWeight(xi0, xi1, xia, xib);
#ifdef DEBUG
                msg = "cellId: " + std::to_string(cellId) + " edgeId: " + std::to_string(edgeIndex) + " weight: " + std::to_string(weight);
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

    loc->Delete();

    return ier;
}

LIBRARY_API
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self, const double data[], double* result) {

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
        std::cout << "adding weight = " << wght << " * edge integral = " << data[cellId*4 + edgeIndex] << " to the flux (" << 
                     *result << " so far) for cellId = " << cellId << " edgeIndex = " << edgeIndex << "\n";
#endif
        *result += wght * data[cellId*4 + edgeIndex];
    }
#ifdef DEBUG
    std::cout << "total flux = " << *result << '\n';
#endif

    return 0;
}
