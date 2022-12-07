#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>

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
int mnt_vectorinterp__getEdgeVectorsFromCellByCellData(VectorInterp_t** self, 
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
int mnt_vectorinterp__getFaceVectorsFromCellByCellData(VectorInterp_t** self, 
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
int mnt_vectorinterp__getEdgeVectorsFromUniqueEdgeData(VectorInterp_t** self, 
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
int mnt_vectorinterp__getFaceVectorsFromUniqueEdgeData(VectorInterp_t** self, 
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
        ier = mnt_vectorinterp__getEdgeVectorsFromCellByCellData(self, data, vectors);
    }
    else {
        ier = mnt_vectorinterp__getEdgeVectorsFromUniqueEdgeData(self, data, vectors);
    }
    return ier;
}

LIBRARY_API
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self, 
                                    const double data[], int placement,
                                    double vectors[]) {
    int ier;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        ier = mnt_vectorinterp__getFaceVectorsFromCellByCellData(self, data, vectors);
    }
    else {
        ier = mnt_vectorinterp__getFaceVectorsFromUniqueEdgeData(self, data, vectors);
    }
    return ier;
}

LIBRARY_API
int mnt_vectorinterp_getVectorsOnEdges(VectorInterp_t** self,
                                            const double data[],
                                            int placement,
                                            double u[], double v[],
                                            int fs) {

    std::string msg;
    int ier;
    int numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;
    }

    const double deg2rad = M_PI/180.;

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) numFailures++;

    std::size_t numEdges;
    ier = mnt_grid_getNumberOfEdges(&(*self)->grid, &numEdges);
    if (ier != 0) numFailures++;

    // initialize to zero
    for (auto i = 0; i < numEdges; ++i) {
        u[i] = 0.0;
        v[i] = 0.0;
    }

    // count the number of cells that are adjacent to each edge
    std::vector<int> numAjdacentCells(numEdges, 0);

    // vertices in radians
    Vec3 v0, v1, v2, v3;
    for (std::size_t cellId = 0; cellId < numCells; ++cellId) {

        // get the vertex coords from the two opposite edges
        ier = mnt_grid_getPoints(&(*self)->grid, cellId, 0, &v0[0], &v1[0]);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getPoints(&(*self)->grid, cellId, 2, &v3[0], &v2[0]);
        if (ier != 0) numFailures++;

        // convert to radians
        v0 *= deg2rad;
        v1 *= deg2rad;
        v2 *= deg2rad;
        v3 *= deg2rad;

        // compute the Cartesian coordinates, assuming a radius of one
        Vec3 xyz0 = cartesianFromRadians(v0);
        Vec3 xyz1 = cartesianFromRadians(v1);
        Vec3 xyz2 = cartesianFromRadians(v2);
        Vec3 xyz3 = cartesianFromRadians(v3);

        // unit normals on each vertex
        Vec3 normal0 = xyz0 / sqrt(dot(xyz0, xyz0));
        Vec3 normal1 = xyz1 / sqrt(dot(xyz1, xyz1));
        Vec3 normal2 = xyz2 / sqrt(dot(xyz2, xyz2));
        Vec3 normal3 = xyz3 / sqrt(dot(xyz3, xyz3));

        // Jacobians on vertices
        const double jac0 = crossDot(xyz1 - xyz0, xyz3 - xyz0, normal0);
        const double jac1 = crossDot(xyz2 - xyz1, xyz0 - xyz1, normal1);
        const double jac2 = crossDot(xyz3 - xyz2, xyz1 - xyz2, normal2);
        const double jac3 = crossDot(xyz0 - xyz3, xyz2 - xyz3, normal3);

        // Jacobians on edge centres are the average between two nodes
        const double j01 = 0.5*(jac0 + jac1);
        const double j12 = 0.5*(jac1 + jac2);
        const double j23 = 0.5*(jac2 + jac3);
        const double j30 = 0.5*(jac3 + jac0);

        // covariant bases at cell centres
        const Vec3 drdXsi = 0.5*(xyz1 - xyz0 + xyz2 - xyz3);
        const Vec3 drdEta = 0.5*(xyz3 - xyz0 + xyz2 - xyz1);

        // edge data index
        std::size_t k0, k1, k2, k3;
        std::size_t edgeId0, edgeId1, edgeId2, edgeId3;
        int edgeSign0, edgeSign1, edgeSign2, edgeSign3;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 0, &edgeId0, &edgeSign0);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 1, &edgeId1, &edgeSign1);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 2, &edgeId2, &edgeSign2);
        if (ier != 0) numFailures++;
        ier = mnt_grid_getEdgeId(&(*self)->grid, cellId, 3, &edgeId3, &edgeSign3);
        if (ier != 0) numFailures++;

        if (placement == MNT_CELL_BY_CELL_DATA) {
            edgeSign0 = 1;
            edgeSign1 = 1;
            edgeSign2 = 1;
            edgeSign3 = 1;
            k0 = cellId*MNT_NUM_EDGES_PER_QUAD + 0;
            k1 = cellId*MNT_NUM_EDGES_PER_QUAD + 1;
            k2 = cellId*MNT_NUM_EDGES_PER_QUAD + 2;
            k3 = cellId*MNT_NUM_EDGES_PER_QUAD + 3;
        }
        else {
            k0 = edgeId0;
            k1 = edgeId1;
            k2 = edgeId2;
            k3 = edgeId3;
        }

        double data0 = data[k0] * edgeSign0;
        double data1 = data[k1] * edgeSign1;
        double data2 = data[k2] * edgeSign2;
        double data3 = data[k3] * edgeSign3;

        Vec3 vec0, vec1, vec2, vec3, adjVecXsi, adjVecEta;
        if (fs == MNT_FUNC_SPACE_W2) {

            // W2

            adjVecXsi = 0.5*(data3/j30 + data1/j12)*drdXsi;
            adjVecEta = 0.5*(data0/j01 + data2/j23)*drdEta;
            // dx cross hat{z} points to the negative direction for edges 0 and 2, hence the 
            // negative sign for edges 0 and 2
            vec0 = - data0*drdEta/j01 + adjVecXsi;
            vec1 = + data1*drdXsi/j12 - adjVecEta;
            vec2 = - data2*drdEta/j23 + adjVecXsi;
            vec3 = + data3*drdXsi/j30 - adjVecEta;
        }
        else {

            // W1

            // normals on mid edge location
            Vec3 normal01 = 0.5*(normal0 + normal1);
            Vec3 normal12 = 0.5*(normal1 + normal2);
            Vec3 normal23 = 0.5*(normal2 + normal3);
            Vec3 normal30 = 0.5*(normal3 + normal0);

            // vectors from adjacent edges, evaluated at the mid edge location
            adjVecXsi = 0.5*( data0*cross(drdEta, normal01)/j01 + data2*cross(drdEta, normal23)/j23 );
            adjVecEta = 0.5*( data3*cross(normal30, drdXsi)/j30 + data1*cross(normal12, drdXsi)/j12 );
            vec0 = data0*cross(drdEta, normal01)/j01 + adjVecEta;
            vec1 = data1*cross(normal12, drdXsi)/j12 + adjVecXsi;
            vec2 = data2*cross(drdEta, normal23)/j23 + adjVecEta;
            vec3 = data3*cross(normal30, drdXsi)/j30 + adjVecXsi;
        }


        double cosLam01 = cos(0.5*(v0[LON_INDEX] + v1[LON_INDEX]));
        double cosLam12 = cos(0.5*(v1[LON_INDEX] + v2[LON_INDEX]));
        double cosLam23 = cos(0.5*(v2[LON_INDEX] + v3[LON_INDEX]));
        double cosLam30 = cos(0.5*(v3[LON_INDEX] + v0[LON_INDEX]));

        double sinLam01 = sin(0.5*(v0[LON_INDEX] + v1[LON_INDEX]));
        double sinLam12 = sin(0.5*(v1[LON_INDEX] + v2[LON_INDEX]));
        double sinLam23 = sin(0.5*(v2[LON_INDEX] + v3[LON_INDEX]));
        double sinLam30 = sin(0.5*(v3[LON_INDEX] + v0[LON_INDEX]));

        double cosThe01 = cos(0.5*(v0[LAT_INDEX] + v1[LAT_INDEX]));
        double cosThe12 = cos(0.5*(v1[LAT_INDEX] + v2[LAT_INDEX]));
        double cosThe23 = cos(0.5*(v2[LAT_INDEX] + v3[LAT_INDEX]));
        double cosThe30 = cos(0.5*(v3[LAT_INDEX] + v0[LAT_INDEX]));

        double sinThe01 = sin(0.5*(v0[LAT_INDEX] + v1[LAT_INDEX]));
        double sinThe12 = sin(0.5*(v1[LAT_INDEX] + v2[LAT_INDEX]));
        double sinThe23 = sin(0.5*(v2[LAT_INDEX] + v3[LAT_INDEX]));
        double sinThe30 = sin(0.5*(v3[LAT_INDEX] + v0[LAT_INDEX]));

        // unit vectors in the lon and lat directions, expressed in Cartesian coords
        Vec3 hatLam01;
        hatLam01[0] = -sinLam01;
        hatLam01[1] = +cosLam01;
        hatLam01[2] = 0.;
        Vec3 hatLam12;
        hatLam12[0] = -sinLam12;
        hatLam12[1] = +cosLam12;
        hatLam12[2] = 0.;
        Vec3 hatLam23;
        hatLam23[0] = -sinLam23;
        hatLam23[1] = +cosLam23;
        hatLam23[2] = 0.;
        Vec3 hatLam30;
        hatLam30[0] = -sinLam30;
        hatLam30[1] = +cosLam30;
        hatLam30[2] = 0.;

        Vec3 hatThe01;
        hatThe01[0] = -sinThe01 * cosLam01;
        hatThe01[1] = -sinThe01 * sinLam01;
        hatThe01[2] = +cosThe01;
        Vec3 hatThe12;
        hatThe12[0] = -sinThe12 * cosLam12;
        hatThe12[1] = -sinThe12 * sinLam12;
        hatThe12[2] = +cosThe12;
        Vec3 hatThe23;
        hatThe23[0] = -sinThe23 * cosLam23;
        hatThe23[1] = -sinThe23 * sinLam23;
        hatThe23[2] = +cosThe23;
        Vec3 hatThe30;
        hatThe30[0] = -sinThe30 * cosLam30;
        hatThe30[1] = -sinThe30 * sinLam30;
        hatThe30[2] = +cosThe30;


        // the data are dimensioned num cells * 4 but the vectors are always dimensions num edges
        u[edgeId0] += dot(vec0, hatLam01);
        v[edgeId0] += dot(vec0, hatThe01);
        numAjdacentCells[edgeId0]++;

        u[edgeId1] += dot(vec1, hatLam12);
        v[edgeId1] += dot(vec1, hatThe12);
        numAjdacentCells[edgeId1]++;

        u[edgeId2] += dot(vec2, hatLam23);
        v[edgeId2] += dot(vec2, hatThe23);
        numAjdacentCells[edgeId2]++;

        u[edgeId3] += dot(vec3, hatLam30);
        v[edgeId3] += dot(vec3, hatThe30);
        numAjdacentCells[edgeId3]++;
    }

    // all the edges that divide two cells have been double counted, now correcting
    for (std::size_t i = 0; i < numEdges; ++i) {

        int count = std::max(1, numAjdacentCells[i]);

        u[i] /= count;
        v[i] /= count;
    }


    return numFailures;

}


