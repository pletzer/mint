#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>

#include "mntLogger.h"
#include <mntExtensiveFieldAdaptor.h>

LIBRARY_API 
int mnt_extensivefieldadaptor_new(ExtensiveFieldAdaptor_t** self) {
     *self = new ExtensiveFieldAdaptor_t();
    (*self)->grid = NULL;
    (*self)->numCells = 0;
    (*self)->numEdges = 0;
    return 0;   
}

LIBRARY_API
int mnt_extensivefieldadaptor_del(ExtensiveFieldAdaptor_t** self) {
    delete *self;
    return 0;
}

LIBRARY_API
int mnt_extensivefieldadaptor_setGrid(ExtensiveFieldAdaptor_t** self, Grid_t* grid) {

    (*self)->grid = grid;

    int ier;

    ier = mnt_grid_getNumberOfCells(&(*self)->grid, &(*self)->numCells);
    if (ier != 0) {
        return ier;
    }

    // ok to have zero if grid was not built from Ugrid
    mnt_grid_getNumberOfEdges(&(*self)->grid, &(*self)->numEdges);

    return 0;    
}

LIBRARY_API
int mnt_extensivefieldadaptor_fromVectorField(ExtensiveFieldAdaptor_t** self,
                                               const double* u, const double* v,
                                              double* data, int placement, int fs) {

    int ier;

    switch (placement) {

        case MNT_CELL_BY_CELL_DATA:

            if (fs == MNT_FUNC_SPACE_W1) {
                ier = mnt_extensivefieldadaptor__fromVectorFieldEdgeCellByCellData(self, u, v, data);
            }
            else {
                ier = mnt_extensivefieldadaptor__fromVectorFieldFaceCellByCellData(self, u, v, data);
            }

            break;

        case MNT_UNIQUE_EDGE_DATA:

            if (fs == MNT_FUNC_SPACE_W1) {
                ier = mnt_extensivefieldadaptor__fromVectorFieldEdgeUniqueIdData(self, u, v, data);
            }
            else {
                ier = mnt_extensivefieldadaptor__fromVectorFieldFaceUniqueIdData(self, u, v, data);
            }

            break;

        default:

            return -1;
    }

    return ier;
}

LIBRARY_API
int mnt_extensivefieldadaptor_toVectorField(ExtensiveFieldAdaptor_t** self, 
                                            const double* edgeData,
                                            const double* faceData,
                                            double* u, double* v,
                                            int placement) {
    int ier;

    switch (placement) {

        case MNT_CELL_BY_CELL_DATA:

            ier = mnt_extensivefieldadaptor__toVectorFieldCellByCellData(self,
                                                    edgeData, faceData, u, v);
            break;

        case MNT_UNIQUE_EDGE_DATA:

            ier = mnt_extensivefieldadaptor__toVectorFieldUniqueIdData(self,
                                                    edgeData, faceData, u, v);
            break;

        default:

            return -1;
    }

    return ier;

}

// private methods
int mnt_extensivefieldadaptor__fromVectorFieldEdgeUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                            const double* u, const double* v,
                                                            double* data) {
    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    if ((*self)->numEdges == 0) {
        msg ="number of edges is zero, the grid was likely not built from UGRID file or data";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;        
    }

    std::size_t edgeId;
    int edgeSign;

    double deg2rad = 0;
    if ((*self)->grid->degrees) {
        deg2rad = M_PI/180.0;
    }

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {

        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getEdgeId(&(*self)->grid, icell, iedge, &edgeId, &edgeSign);
            if (ier != 0) numFailures++;
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double x0 = point0[LON_INDEX];
            double x1 = point1[LON_INDEX];
            double dx = x1 - x0;
            double y0 = point0[LAT_INDEX];
            double y1 = point1[LAT_INDEX];
            double dy = y1 - y0;

            // length on the surface of the sphere (1 if no spherical)
            double cosTheta = cos( deg2rad * 0.5 * (y1 + y0) );
            double cosTheDx = cosTheta * dx;

            // edge integral
            data[edgeId] = edgeSign * ( u[edgeId] * cosTheDx + v[edgeId] * dy );
        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__fromVectorFieldFaceUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                            const double* u, const double* v,
                                                            double* data) {
    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    if ((*self)->numEdges == 0) {
        msg ="number of edges is zero, the grid was likely not built from UGRID file or data";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;        
    }

    std::size_t edgeId;
    int edgeSign;

    double deg2rad = 0;
    if ((*self)->grid->degrees) {
        deg2rad = M_PI/180.0;
    }

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {

        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getEdgeId(&(*self)->grid, icell, iedge, &edgeId, &edgeSign);
            if (ier != 0) numFailures++;
        
            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double x0 = point0[LON_INDEX];
            double x1 = point1[LON_INDEX];
            double dx = x1 - x0;
            double y0 = point0[LAT_INDEX];
            double y1 = point1[LAT_INDEX];
            double dy = y1 - y0;

            // length on the surface of the sphere (1 if no spherical)
            double cosTheta = cos( deg2rad * 0.5 * (y1 + y0) );
            double cosTheDx = cosTheta * dx;

            // flux integral

            // our convention is to have the extensive fluxes pointing in the positive, logical direction.
            // This requires flipping the sign for iedge = 0 and 2.
            double sign = -2*((iedge + 1) % 2) + 1;

            data[edgeId] = edgeSign*sign*(u[edgeId]*dy - v[edgeId]*cosTheDx);
        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__fromVectorFieldEdgeCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double* u, const double* v,
                                                            double* data) {
    std::string msg = "NOT IMPLEMENTED";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return -1;
}

int mnt_extensivefieldadaptor__fromVectorFieldFaceCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double* u, const double* v,
                                                            double* data) {
    std::string msg = "NOT IMPLEMENTED";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return -1;
}

int mnt_extensivefieldadaptor__toVectorFieldUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                           const double* edgeData,
                                                           const double* faceData,
                                                           double* u, double* v) {
    std::string msg;
    int ier, numFailures = 0;

    if (!(*self)->grid) {
        msg ="must set the grid before calling this";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -1;
    }

    if ((*self)->numEdges == 0) {
        msg ="number of edges is zero, the grid was likely not built from UGRID file or data";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return -2;        
    }

    std::size_t edgeId;
    int edgeSign;

    double deg2rad = 0;
    if ((*self)->grid->degrees) {
        deg2rad = M_PI/180.0;
    }

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {

        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getEdgeId(&(*self)->grid, icell, iedge, &edgeId, &edgeSign);
            if (ier != 0) numFailures++;

            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double x0 = point0[LON_INDEX];
            double x1 = point1[LON_INDEX];
            double dx = x1 - x0;
            double y0 = point0[LAT_INDEX];
            double y1 = point1[LAT_INDEX];
            double dy = y1 - y0;

            // length on the surface of the sphere (1 if no spherical)
            double cosTheta = cos( deg2rad * 0.5 * (y1 + y0) );
            double cosTheDx = cosTheta * dx;

            // our convention is to have the extensive fluxes pointing in the positive, logical direction.
            // This requires flipping the sign for iedge = 0 and 2.
            double sign = -2*((iedge + 1) % 2) + 1;

            double len_sq = cosTheDx * cosTheDx + dy * dy + 1.e-10; // allow for a small epsilon at the poles

            u[edgeId] = edgeSign * (cosTheDx * edgeData[edgeId] + sign * dy * faceData[edgeId] ) / len_sq;
            v[edgeId] = edgeSign * (dy * edgeData[edgeId] - sign * cosTheDx * faceData[edgeId] ) / len_sq;
        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__toVectorFieldCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                           const double* edgeData,
                                                           const double* faceData,
                                                           double* u, double* v) {
    std::string msg = "NOT IMPLEMENTED";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return -1;
}
