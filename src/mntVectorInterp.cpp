#include <mntVectorInterp.h>
#include <vtkGenericCell.h>


extern "C" 
int mnt_vectorinterp_new(VectorInterp_t** self) {
    *self = new VectorInterp_t();
    (*self)->grid = NULL;
    (*self)->locator = NULL;
    (*self)->ownsLocator = false;
    (*self)->calledFindPoints = false;
    return 0;
}

extern "C"
int mnt_vectorinterp_del(VectorInterp_t** self) {

    if ((*self)->ownsLocator) {
        (*self)->locator->Delete();
    }

    delete *self;

    return 0;
}

extern "C"
int mnt_vectorinterp_setGrid(VectorInterp_t** self, Grid_t* grid) {
    (*self)->grid = grid;
    return 0;
}

extern "C"
int mnt_vectorinterp_setLocator(VectorInterp_t** self, vmtCellLocator* locator) {
    (*self)->locator = locator;
    return 0;
}

extern "C"
int mnt_vectorinterp_buildLocator(VectorInterp_t** self, int numCellsPerBucket, double periodX) {

    if (!(*self)->grid) {
        std::cerr << "ERROR: must call setGrid before invoking buildLocator\n";
        return -1;
    }

    (*self)->ownsLocator = true;
    (*self)->locator = vmtCellLocator::New();
    (*self)->locator->SetDataSet((*self)->grid->grid);
    (*self)->locator->SetNumberOfCellsPerBucket(numCellsPerBucket);
    (*self)->locator->BuildLocator();
    (*self)->locator->setPeriodicityLengthX(periodX);

    return 0;
}

extern "C"
int mnt_vectorinterp_findPoints(VectorInterp_t** self, std::size_t numPoints, 
                                    const double targetPoints[], double tol2) {

    if (!(*self)->locator) {
        std::cerr << "ERROR: must call either setLocator or buildLocator before invoking findPoints\n";
        return -1;
    }

    vtkGenericCell *cell = NULL;
    double weights[8];
    int numFailures = 0;
    double pcoords[3];

    (*self)->cellIds.resize(numPoints);
    (*self)->pcoords.resize(numPoints);

    for (std::size_t i = 0; i < numPoints; ++i) {

        // find the cell
        (*self)->cellIds[i] = (*self)->locator->FindCell(&targetPoints[3*i], tol2,
                                                         cell, pcoords, weights);
        if ((*self)->cellIds[i] < 0) {
            // outside the domain?
            std::cerr << "Warning: target point " << 
                          targetPoints[3*i] << ',' << targetPoints[3*i + 1] << ',' << targetPoints[3*i + 2] << 
                          " is outside of the domain\n";
            numFailures++;
            continue;
        }

        // store the parametric coordinates of the target point
        for (std::size_t j = 0; j < 3; ++j) {
            (*self)->pcoords[i][j] = pcoords[j];
        }

    }

    (*self)->calledFindPoints = true;

    return numFailures;
}

extern "C"
int mnt_vectorinterp_getEdgeVectors(VectorInterp_t** self, 
                                    const double data[], double vectors[]) {

    if (!(*self)->calledFindPoints) {
        std::cerr << "ERROR: must call findPoints before calling getEdgeVectors\n";
        return -1;
    }
    if (!(*self)->grid) {
        std::cerr << "ERROR: must set the grid before calling getEdgeVectors\n";
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        std::cerr << "Warning: there is no target point in the domain\n";
        return 2;
    }

    int numFailures = 0;
    Vec3 drdXsi, drdEta, gradXsi, gradEta;
    double jac;

    for (std::size_t iTargetId = 0;
                     iTargetId < (*self)->cellIds.size(); 
                     ++iTargetId) {

        vtkIdType cellId = (*self)->cellIds[iTargetId];

        if (cellId < 0) {
            numFailures++;
            // next target
            continue;
        }

        mnt_vectorinterp__getTangentVectors(self, iTargetId, drdXsi, drdEta, jac);

        // contravariant bases
        gradXsi[0] = + drdEta[1]/jac;
        gradXsi[1] = - drdEta[0]/jac;
        gradXsi[2] = 0.0;

        gradEta[0] = - drdXsi[1]/jac;
        gradEta[1] = + drdXsi[0]/jac;
        gradEta[2] = 0.0;

        // parametric coordinates of the target point 
        double xsi = (*self)->pcoords[iTargetId][0];
        double eta = (*self)->pcoords[iTargetId][1];
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;

        // interpolate
        double data0 = data[cellId*4 + 0];
        double data1 = data[cellId*4 + 1];
        double data2 = data[cellId*4 + 2];
        double data3 = data[cellId*4 + 3];
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            vectors[iTargetId*3 + j] = (data0*ate + data2*eta)*gradXsi[j] + 
                                       (data3*isx + data1*xsi)*gradEta[j];
        }
    }

    return numFailures;
}

extern "C"
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self, 
                                    const double data[], double vectors[]) {

    if (!(*self)->calledFindPoints) {
        std::cerr << "ERROR: must call findPoints before calling getFaceVectors\n";
        return -1;
    }
    if (!(*self)->grid) {
        std::cerr << "ERROR: must set the grid before calling getFaceVectors\n";
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        std::cerr << "Warning: there is no target point in the domain\n";
        return 2;
    }

    int numFailures = 0;
    Vec3 drdXsi, drdEta, gradXsi, gradEta;
    double jac;

    for (std::size_t iTargetId = 0;
                     iTargetId < (*self)->cellIds.size(); 
                     ++iTargetId) {

        vtkIdType cellId = (*self)->cellIds[iTargetId];

        if (cellId < 0) {
            numFailures++;
            // next target
            continue;
        }

        mnt_vectorinterp__getTangentVectors(self, iTargetId, drdXsi, drdEta, jac);

        // parametric coordinates of the target point 
        double xsi = (*self)->pcoords[iTargetId][0];
        double eta = (*self)->pcoords[iTargetId][1];
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;

        // interpolate
        double data0 = data[cellId*4 + 0];
        double data1 = data[cellId*4 + 1];
        double data2 = data[cellId*4 + 2];
        double data3 = data[cellId*4 + 3];
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            // the negative sign is because d r/ d eta points in the opposite
            /// direction to the surface
            vectors[iTargetId*3 + j] = (data3*isx + data1*xsi)*drdXsi[j]/jac 
                                     - (data0*ate + data2*eta)*drdEta[j]/jac;
        }
    }

    return numFailures;
}

