#include <mntExtensiveFieldConverter.h>

LIBRARY_API 
int mnt_extensivefieldconverter_new(ExtensiveFieldConverter_t** self) {
    *self = new ExtensiveFieldConverter_t();
    (*self)->grid = NULL;
    return 0;
}

LIBRARY_API
int mnt_extensivefieldconverter_del(ExtensiveFieldConverter_t** self) {
    delete *self;
    return 0;
}

LIBRARY_API
int mnt_extensivefieldconverter_setGrid(ExtensiveFieldConverter_t** self, Grid_t* grid) {
    (*self)->grid = grid;
    return 0;
}

LIBRARY_API
int mnt_extensivefieldconverter_getEdgeDataFromCellByCellVectors(ExtensiveFieldConverter_t** self, 
                                                      const double vx[], const double vy[],
                                                      double data[]) {

    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) numFailures++;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) numCells; ++icell) {
        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {
            std::size_t k = icell * MNT_NUM_EDGES_PER_QUAD + iedge;
            double u = vx[k];
            double v = vy[k];
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double dx = point1[LON_INDEX] - point0[LON_INDEX];
            double ymid = 0.5*(point1[LAT_INDEX] + point0[LAT_INDEX]);
            double dy = point1[LAT_INDEX] - point0[LAT_INDEX];

            // line integral 
            data[k] = u*dx + v*dy;
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_extensivefieldconverter_getEdgeDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self, 
                                                      const double vx[], const double vy[],
                                                      double data[]) {

    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    std::size_t numEdges = 0;
    ier = mnt_grid_getNumberOfEdges(&(*self)->grid, &numEdges);
    if (ier !=0 || numEdges == 0) {
        msg ="number of edges is zero, the grid was likely not built from a UGRID file";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;        
    }

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) numFailures++;

    std::size_t edgeId;
    int edgeSign;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) numCells; ++icell) {
        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getEdgeId(&(*self)->grid, icell, iedge, &edgeId, &edgeSign);
            if (ier != 0) numFailures++;

            // unique edges
            double u = vx[edgeId];
            double v = vy[edgeId];
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double dx = point1[LON_INDEX] - point0[LON_INDEX];
            double dy = point1[LAT_INDEX] - point0[LAT_INDEX];

            // data is cell by cell even though vx and vy are on unique edges
            std::size_t k = icell * MNT_NUM_EDGES_PER_QUAD + iedge;
            data[k] = u*dx + v*dy;
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_extensivefieldconverter_getFaceDataFromCellByCellVectors(ExtensiveFieldConverter_t** self, 
                                                      const double vx[], const double vy[],
                                                      double data[]) {

    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) numFailures++;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) numCells; ++icell) {
        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {
            std::size_t k = icell * MNT_NUM_EDGES_PER_QUAD + iedge;
            double u = vx[k];
            double v = vy[k];
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double dx = point1[LON_INDEX] - point0[LON_INDEX];
            double dy = point1[LAT_INDEX] - point0[LAT_INDEX];

            // line integral 
            double sign = 2*((iedge + 1) % 2) - 1;
            data[k] = -sign*(u*dy - v*dx);
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_extensivefieldconverter_getFaceDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self, 
                                                      const double vx[], const double vy[],
                                                      double data[]) {

    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    std::size_t numEdges = 0;
    ier = mnt_grid_getNumberOfEdges(&(*self)->grid, &numEdges);
    if (ier !=0 || numEdges == 0) {
        msg ="number of edges is zero, the grid was likely not built from a UGRID file";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;        
    }

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &numCells);
    if (ier != 0) {
        numFailures++;
    }

    std::size_t edgeId;
    int edgeSign;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) numCells; ++icell) {
        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getEdgeId(&(*self)->grid, icell, iedge, &edgeId, &edgeSign);
            if (ier != 0) numFailures++;

            // unique edges
            double u = vx[edgeId];
            double v = vy[edgeId];
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double dx = point1[LON_INDEX] - point0[LON_INDEX];
            double dy = point1[LAT_INDEX] - point0[LAT_INDEX];

            // data is cell by cell even though vx and vy are on unique edges
            std::size_t k = icell * MNT_NUM_EDGES_PER_QUAD + iedge;
            double sign = 2*((iedge + 1) % 2) - 1;
            data[k] = -sign*(u*dy - v*dx);
        }
    }

    return numFailures;
}

LIBRARY_API
int mnt_extensivefieldconverter_getEdgeData(ExtensiveFieldConverter_t** self,
                                            const double vx[], const double vy[], int placement,
                                            double data[]) {
    int ier;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        ier = mnt_extensivefieldconverter_getEdgeDataFromCellByCellVectors(self, vx, vy, data);
    }
    else {
        ier = mnt_extensivefieldconverter_getEdgeDataFromUniqueEdgeVectors(self, vx, vy, data);
    }
    return ier;
}

LIBRARY_API
int mnt_extensivefieldconverter_getFaceData(ExtensiveFieldConverter_t** self,
                                            const double vx[], const double vy[], int placement,
                                            double data[]) {
    int ier;
    if (placement == MNT_CELL_BY_CELL_DATA) {
        ier = mnt_extensivefieldconverter_getFaceDataFromCellByCellVectors(self, vx, vy, data);
    }
    else {
        ier = mnt_extensivefieldconverter_getFaceDataFromUniqueEdgeVectors(self, vx, vy, data);
    }
    return ier;
}

