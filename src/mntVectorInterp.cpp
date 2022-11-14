#include <mntVectorInterp.h>
#include <vtkGenericCell.h>


LIBRARY_API 
int mnt_vectorinterp_new(VectorInterp_t** self) {
    *self = new VectorInterp_t();
    (*self)->grid = NULL;
    (*self)->locator = NULL;
    (*self)->ownsLocator = false;
    (*self)->calledFindPoints = false;
    return 0;
}

LIBRARY_API
int mnt_vectorinterp_del(VectorInterp_t** self) {

    if ((*self)->ownsLocator) {
        (*self)->locator->Delete();
    }

    delete *self;

    return 0;
}

LIBRARY_API
int mnt_vectorinterp_setGrid(VectorInterp_t** self, Grid_t* grid) {
    (*self)->grid = grid;
    return 0;
}

LIBRARY_API
int mnt_vectorinterp_setLocator(VectorInterp_t** self, vmtCellLocator* locator) {
    (*self)->locator = locator;
    return 0;
}

LIBRARY_API
int mnt_vectorinterp_buildLocator(VectorInterp_t** self, int numCellsPerBucket, double periodX, int enableFolding) {

    if (!(*self)->grid) {
        std::string msg ="must call setGrid before invoking buildLocator";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    (*self)->ownsLocator = true;
    (*self)->locator = vmtCellLocator::New();
    (*self)->locator->SetDataSet((*self)->grid->grid);
    (*self)->locator->SetNumberOfCellsPerBucket(numCellsPerBucket);
    (*self)->locator->setPeriodicityLengthX(periodX);
    if (enableFolding == 1) {
        (*self)->locator->enableFolding();
    }
    (*self)->locator->BuildLocator();

    return 0;
}

LIBRARY_API
int mnt_vectorinterp_findPoints(VectorInterp_t** self, std::size_t numPoints, 
                                    const double targetPoints[], double tol2) {

    std::string msg;

    if (!(*self)->locator) {
        msg ="must call either setLocator or buildLocator before invoking findPoints";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
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
            msg ="target point " + 
                 std::to_string(targetPoints[3*i    ]) + ',' + 
                 std::to_string(targetPoints[3*i + 1]) + ',' + 
                 std::to_string(targetPoints[3*i + 2]) + 
                 " is outside of the domain";
            mntlog::warn(__FILE__, __func__, __LINE__, msg);
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

LIBRARY_API
int mnt_vectorinterp_getEdgeVectorsFromCellByCellData(VectorInterp_t** self, 
                                                      const double data[],
                                                      double vectors[]) {

    std::string msg;

    if (!(*self)->calledFindPoints) {
        msg ="must call findPoints before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        msg ="there is no target point in the domain";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
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
        double data0 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 0];
        double data1 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 1];
        double data2 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 2];
        double data3 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 3];
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            vectors[iTargetId*3 + j] = (data0*ate + data2*eta)*gradXsi[j] + 
                                       (data3*isx + data1*xsi)*gradEta[j];
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_vectorinterp_getFaceVectorsFromCellByCellData(VectorInterp_t** self, 
                                                      const double data[],
                                                      double vectors[]) {

    std::string msg;

    if (!(*self)->calledFindPoints) {
        msg ="must call findPoints before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        msg ="there is no target point in the domain";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
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
        double data0 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 0];
        double data1 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 1];
        double data2 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 2];
        double data3 = data[cellId*MNT_NUM_EDGES_PER_QUAD + 3];
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

LIBRARY_API
int mnt_vectorinterp_getEdgeVectorsFromUniqueEdgeData(VectorInterp_t** self, 
                                                      const double data[],
                                                      double vectors[]) {
    
    std::string msg;
    std::size_t edgeId0, edgeId1, edgeId2, edgeId3;
    int edgeSign0, edgeSign1, edgeSign2, edgeSign3;

    if (!(*self)->calledFindPoints) {
        msg ="must call findPoints before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        msg ="there is no target point in the domain";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
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

        mnt_grid_getEdgeId(&(*self)->grid, cellId, 0, &edgeId0, &edgeSign0);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 1, &edgeId1, &edgeSign1);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 2, &edgeId2, &edgeSign2);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 3, &edgeId3, &edgeSign3);

        // interpolate
        double data0 = data[edgeId0] * edgeSign0;
        double data1 = data[edgeId1] * edgeSign1;
        double data2 = data[edgeId2] * edgeSign2;
        double data3 = data[edgeId3] * edgeSign3;
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            vectors[iTargetId*3 + j] = (data0*ate + data2*eta)*gradXsi[j] + 
                                       (data3*isx + data1*xsi)*gradEta[j];
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_vectorinterp_getFaceVectorsFromUniqueEdgeData(VectorInterp_t** self, 
                                                      const double data[],
                                                      double vectors[]) {

    std::string msg;
    std::size_t edgeId0, edgeId1, edgeId2, edgeId3;
    int edgeSign0, edgeSign1, edgeSign2, edgeSign3;

    if (!(*self)->calledFindPoints) {
        msg ="must call findPoints before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }
    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    if ((*self)->pcoords.size() == 0) {
        msg ="there is no target point in the domain";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
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

        mnt_grid_getEdgeId(&(*self)->grid, cellId, 0, &edgeId0, &edgeSign0);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 1, &edgeId1, &edgeSign1);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 2, &edgeId2, &edgeSign2);
        mnt_grid_getEdgeId(&(*self)->grid, cellId, 3, &edgeId3, &edgeSign3);

        // interpolate
        double data0 = data[edgeId0] * edgeSign0;
        double data1 = data[edgeId1] * edgeSign1;
        double data2 = data[edgeId2] * edgeSign2;
        double data3 = data[edgeId3] * edgeSign3;
        for (std::size_t j = 0; j < 3; ++j) {
            // fill in the vector
            // the negative sign is because d r/ d eta points in the opposite
            // direction to the surface
            vectors[iTargetId*3 + j] = (data3*isx + data1*xsi)*drdXsi[j]/jac 
                                     - (data0*ate + data2*eta)*drdEta[j]/jac;
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_vectorinterp_getEdgeVectors(VectorInterp_t** self, 
                                    const double data[], int placement,
                                    double vectors[]) {
    int ier;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        ier = mnt_vectorinterp_getEdgeVectorsFromCellByCellData(self, data, vectors);
    }
    else {
        ier = mnt_vectorinterp_getEdgeVectorsFromUniqueEdgeData(self, data, vectors);
    }
    return ier;
}

LIBRARY_API
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self, 
                                    const double data[], int placement,
                                    double vectors[]) {
    int ier;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        ier = mnt_vectorinterp_getFaceVectorsFromCellByCellData(self, data, vectors);
    }
    else {
        ier = mnt_vectorinterp_getFaceVectorsFromUniqueEdgeData(self, data, vectors);
    }
    return ier;
}


LIBRARY_API
int mnt_vectorinterp_getFaceVectorsFromUniqueEdgeDataOnEdges(VectorInterp_t** self,
                                                          const double data[],
                                                          double u[], double v[]) {

    std::string msg;
    int ier;
    int numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    Vec3 drdXsi, drdEta;
    std::size_t edgeId0, edgeId1, edgeId2, edgeId3;
    int edgeSign0, edgeSign1, edgeSign2, edgeSign3;

    std::size_t numCells;
    ier =  mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) numFailures++;

    std::size_t numEdges;
    ier =  mnt_grid_getNumberOfEdges(&(*self)->grid, &numEdges);
    if (ier != 0) numFailures++;

    // count the number of cells that are adjacent to each edge
    std::vector<int> numAjdacentCells(numEdges, 0);

    // vertex positions
    Vec3 v0, v1, v2, v3;

    // initialize the vector field to zero
    for (auto i = 0; i < numEdges; ++i) {
        u[i] = 0.0;
        v[i] = 0.0;
    }


    for (std::size_t cellId = 0; cellId < numCells; ++cellId) {

        // get the vertex coords
        ier = mnt_grid_getPoints(&(*self)->grid, cellId, 0, &v0[0], &v1[0]);
        if (ier != 0) numFailures++;

        ier = mnt_grid_getPoints(&(*self)->grid, cellId, 2, &v3[0], &v2[0]);
        if (ier != 0) numFailures++;

        Vec3 a = v1 - v0;
        Vec3 b = v2 - v1;
        Vec3 c = v2 - v3;
        Vec3 d = v3 - v0;

        // Jacobians attached to each vertex (can be zero if points are degenerate)
        double a013 = crossDotZHat(a, d);
        double a120 = crossDotZHat(a, b);
        double a231 = crossDotZHat(c, b);
        double a302 = crossDotZHat(c, d);

        // Jacobian for this quad, should be a strictly positive quantity if nodes are
        // ordered correctly
        double jac = 0.25*(a013 + a120 + a231 + a302);
        if (jac <= 0) {
            std::stringstream msg;
            msg << "bad cell " << cellId << " vertices: " <<
                            v0 << ";" << v1 << ";" << v2  << ";" << v3; 
            mntlog::warn(__FILE__, __func__, __LINE__, msg.str());
            ier = 1;
        }


        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 0, &edgeId0, &edgeSign0);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 1, &edgeId1, &edgeSign1);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 2, &edgeId2, &edgeSign2);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 3, &edgeId3, &edgeSign3);
        if (ier != 0) numFailures++;

        // cotangent vectors obtained by finite differencing and linearly interpolating
        drdXsi = 0.5*(a + c);
        drdEta = 0.5*(d + b);


        // interpolate
        double data0 = data[edgeId0] * edgeSign0;
        double data1 = data[edgeId1] * edgeSign1;
        double data2 = data[edgeId2] * edgeSign2;
        double data3 = data[edgeId3] * edgeSign3;

        // edge 0
        double xsi = 0.5;
        double eta = 0.0;
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;
        u[edgeId0] += (data3*isx + data1*xsi)*drdXsi[0]/jac 
                    - (data0*ate + data2*eta)*drdEta[0]/jac;
        v[edgeId0] += (data3*isx + data1*xsi)*drdXsi[1]/jac 
                    - (data0*ate + data2*eta)*drdEta[1]/jac;
        numAjdacentCells[edgeId0]++;

        // edge 1
        xsi = 1.0;
        eta = 0.5;
        isx = 1.0 - xsi;
        ate = 1.0 - eta;
        u[edgeId1] += (data3*isx + data1*xsi)*drdXsi[0]/jac 
                    - (data0*ate + data2*eta)*drdEta[0]/jac;
        v[edgeId1] += (data3*isx + data1*xsi)*drdXsi[1]/jac 
                    - (data0*ate + data2*eta)*drdEta[1]/jac;
        numAjdacentCells[edgeId1]++;

        // edge 2
        xsi = 0.5;
        eta = 1.0;
        isx = 1.0 - xsi;
        ate = 1.0 - eta;
        u[edgeId2] += (data3*isx + data1*xsi)*drdXsi[0]/jac 
                    - (data0*ate + data2*eta)*drdEta[0]/jac;
        v[edgeId2] += (data3*isx + data1*xsi)*drdXsi[1]/jac 
                    - (data0*ate + data2*eta)*drdEta[1]/jac;
        numAjdacentCells[edgeId2]++;

        // edge 3
        xsi = 0.0;
        eta = 0.5;
        isx = 1.0 - xsi;
        ate = 1.0 - eta;
        u[edgeId3] += (data3*isx + data1*xsi)*drdXsi[0]/jac 
                    - (data0*ate + data2*eta)*drdEta[0]/jac;
        v[edgeId3] += (data3*isx + data1*xsi)*drdXsi[1]/jac 
                    - (data0*ate + data2*eta)*drdEta[1]/jac;
        numAjdacentCells[edgeId3]++;


    }

    // all the edges that divide two cells have been double counted, now correcting
    for (auto i = 0; i < numEdges; ++i) {
        int count = std::max(1, numAjdacentCells[i]);
        u[i] /= count;
        v[i] /= count;
    }

    return numFailures;
}


