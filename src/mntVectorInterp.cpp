#include <mntVectorInterp.h>
#include <vtkGenericCell.h>

double crossDotZHat(const Vec3& a, const Vec3& b) {
    return a[0]*b[1] - a[1]*b[0];
}

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
        std::cerr << "ERROR: must call findPoints before getting the vector\n";
        return -1;
    }
    if (!(*self)->grid) {
        std::cerr << "ERROR: must set the grid\n";
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        std::cerr << "Warning: there is no target point in the domain\n";
        return 2;
    }

    int numFailures = 0;
    Vec3 v0, v1, v2, v3, gradXsi, gradEta;

    for (std::size_t iTargetId = 0;
                     iTargetId < (*self)->cellIds.size(); 
                     ++iTargetId) {

        vtkIdType cellId = (*self)->cellIds[iTargetId];

        if (cellId < 0) {
            numFailures++;
            // next target
            continue;
        }

        // parametric coordinates of the target point 
        double xsi = (*self)->pcoords[iTargetId][0];
        double eta = (*self)->pcoords[iTargetId][1];
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;

        // get the cell vertices, this should never fail 
        mnt_grid_getPoints(&(*self)->grid, cellId, 0, &v0[0], &v1[0]);
        mnt_grid_getPoints(&(*self)->grid, cellId, 2, &v3[0], &v2[0]);

        // Jacobians attached to each vertex (can be zero if points are degenerate)
        double a013 = crossDotZHat(v1 - v0, v3 - v0);
        double a120 = crossDotZHat(v2 - v1, v0 - v1);
        double a231 = crossDotZHat(v3 - v2, v1 - v2);
        double a302 = crossDotZHat(v0 - v3, v2 - v3);

        // Jacobian for this quad, should be a strictly positive quantity if nodes are
        // ordered correctly
        double jac = 0.25*(a013 + a120 + a231 + a302);

        // cotangent vectors obtained by finite differencing and linearly interpolating
        // in the other direction
        Vec3 drdXsi = ate*(v1 - v0) + eta*(v2 - v3);
        Vec3 drdEta = isx*(v3 - v0) + xsi*(v2 - v1);

        // contravariant bases
        gradXsi[0] = + drdEta[1]/jac;
        gradXsi[1] = - drdEta[0]/jac;
        gradXsi[2] = 0.0;

        gradEta[0] = - drdXsi[1]/jac;
        gradEta[1] = + drdXsi[0]/jac;
        gradEta[2] = 0.0;

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

    return 0;
}

extern "C"
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self, 
                                    const double data[], double vectors[]) {

    if (!(*self)->calledFindPoints) {
        std::cerr << "ERROR: must call findPoints before getting the vector\n";
        return -1;
    }
    if (!(*self)->grid) {
        std::cerr << "ERROR: must set the grid\n";
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        std::cerr << "Warning: there is no target point in the domain\n";
        return 2;
    }

    int numFailures = 0;
    Vec3 v0, v1, v2, v3, gradXsi, gradEta;

    for (std::size_t iTargetId = 0;
                     iTargetId < (*self)->cellIds.size(); 
                     ++iTargetId) {

        vtkIdType cellId = (*self)->cellIds[iTargetId];

        if (cellId < 0) {
            numFailures++;
            // next target
            continue;
        }

        // parametric coordinates of the target point 
        double xsi = (*self)->pcoords[iTargetId][0];
        double eta = (*self)->pcoords[iTargetId][1];
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;

        // get the cell vertices, this should never fail 
        mnt_grid_getPoints(&(*self)->grid, cellId, 0, &v0[0], &v1[0]);
        mnt_grid_getPoints(&(*self)->grid, cellId, 2, &v3[0], &v2[0]);

        // Jacobians attached to each vertex (can be zero if points are degenerate)
        double a013 = crossDotZHat(v1 - v0, v3 - v0);
        double a120 = crossDotZHat(v2 - v1, v0 - v1);
        double a231 = crossDotZHat(v3 - v2, v1 - v2);
        double a302 = crossDotZHat(v0 - v3, v2 - v3);

        // Jacobian for this quad, should be a strictly positive quantity if nodes are
        // ordered correctly
        double jac = 0.25*(a013 + a120 + a231 + a302);

        // cotangent vectors obtained by finite differencing and linearly interpolating
        // in the other direction
        Vec3 drdXsi = ate*(v1 - v0) + eta*(v2 - v3);
        Vec3 drdEta = isx*(v3 - v0) + xsi*(v2 - v1);

        // interpolate
        double data0 = data[cellId*4 + 0];
        double data1 = data[cellId*4 + 1];
        double data2 = data[cellId*4 + 2];
        double data3 = data[cellId*4 + 3];
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            vectors[iTargetId*3 + j] = (data0*ate + data2*eta)*drdEta[j]/jac + 
                                       (data3*isx + data1*xsi)*drdXsi[j]/jac;
        }
    }

    return 0;
}

