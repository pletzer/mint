#include <mntPolylineIntegral.h>
#include <mntPolysegmentIter.h>
#include <iostream>

#define DEBUG

extern "C"
int mnt_polylineintegral_new(PolylineIntegral_t** self) {
    
    *self = new PolylineIntegral_t();
    (*self)->grid = NULL;

    return 0;
}

extern "C"
int mnt_polylineintegral_del(PolylineIntegral_t** self) {
    int ier = 0;

    // destroy everything here

    delete *self;
    return ier;
}

extern "C"
int mnt_polylineintegral_setGrid(PolylineIntegral_t** self, Grid_t* grid) {
    (*self)->grid = grid;
    return 0;
}

extern "C"
int mnt_polylineintegral_setPolyline(PolylineIntegral_t** self, int npoints, const double xyz[]) {
    (*self)->lineXYZ.resize(3 * npoints);
    // copy the line
    for (auto i = 0; i < npoints*3; ++i) {
        (*self)->lineXYZ[i] = xyz[i]; // copy
    }
    return 0;
}

extern "C"
int mnt_polylineintegral_build(PolylineIntegral_t** self) {

    if ((*self)->lineXYZ.size() == 0) {
        std::cerr << "ERROR: need to call setPolyline before calling build\n";
        return 1;
    }

    if (!(*self)->grid) {
        std::cerr << "ERROR: need to call setGrid before calling build\n";
        return 2;
    }

    vtkUnstructuredGrid* vgrid = (*self)->grid->grid;

    // build the cell locator
    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(vgrid);
    loc->BuildLocator();

    // number of points/vertices (3D)
    size_t npts = (*self)->lineXYZ.size() / 3;

    // iterate over the Polyline segments
    for (size_t ip0 = 0; ip0 < npts - 1; ++ip0) {

        // get the starting point
        size_t k0 = 3*ip0;
        double p0[] = {(*self)->lineXYZ[k0 + 0], (*self)->lineXYZ[k0 + 1], (*self)->lineXYZ[k0 + 2]};

        // get the end point
        size_t k1 = k0 + 3;        
        double p1[] = {(*self)->lineXYZ[k1 + 0], (*self)->lineXYZ[k1 + 1], (*self)->lineXYZ[k1 + 2]};

        // build the polysegment iterator. This breaks segment p0 -> p1 into 
        // smaller segments, each being entirely contained within grid cells
        PolysegmentIter polyseg(vgrid, loc, p0, p1);

        // number of sub-segments
        size_t numSubSegs = polyseg.getNumberOfSegments();
        polyseg.reset();

        // iterate over the sub-segments 
        for (size_t i = 0; i < numSubSegs; ++i) {

            // get the cell Id to which this sub-segment belongs to
            vtkIdType cellId = polyseg.getCellId();

            // get the cell parametric start/end values for the sub-segment
            const Vec3& xia = polyseg.getBegCellParamCoord();
            const Vec3& xib = polyseg.getEndCellParamCoord();

            // some sub-segments belong to two or more cells (e.g. if the sub-segment
            // runs along a cell edge)
            double coeff = polyseg.getCoefficient();

            // 2D because our cells are quads
            std::vector<double> dxi({xib[0] - xia[0], xib[1] - xia[1]});
            std::vector<double> xiMid({0.5*(xia[0] + xib[0]), 0.5*(xia[1] + xib[1])});

            // compute the 4 weight contributions from each cell edge
            double ws[] = {+ dxi[0] * (1.0 - xiMid[1]) * coeff,
                           + dxi[1] * (0.0 + xiMid[0]) * coeff,
                           - dxi[0] * (0.0 + xiMid[1]) * coeff,
                           - dxi[1] * (1.0 - xiMid[0]) * coeff};

            // store the weights for this sub-segment
            std::pair<vtkIdType, int> cIdE;
            for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex) {
                cIdE = std::pair<vtkIdType, int>(cellId, edgeIndex);
#ifdef DEBUG
                std::cout << "cellId: " << cellId << " edgeIndex: " << edgeIndex << " weight: " << ws[edgeIndex] << '\n';
#endif
                (*self)->weights.insert(std::pair< std::pair<vtkIdType, int>, double>(cIdE, ws[edgeIndex]));
            }            
        }
    }

    loc->Delete();

    return 0;
}

extern "C"
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self, const double data[], double* result) {

    *result = 0;
    for (auto it = (*self)->weights.cbegin(); it != (*self)->weights.cend(); ++it) {

        // get the cell Id
        vtkIdType cellId = it->first.first;

        // get the edge index
        int edgeIndex = it->first.second;

        // get the weight
        double wght = it->second;

        // add the contibution
        *result += wght * data[cellId*4 + edgeIndex];
    }

    return 0;
}
