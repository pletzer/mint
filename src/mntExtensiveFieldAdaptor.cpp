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
                                               const double u[], const double v[],
                                               double data[], int placement, int fs) {

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
                                            const double edgeData[],
                                            const double faceData[],
                                            double u[], double v[],
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

            // length on the surface of the sphere (1 if not spherical)
            double cosTheta = cos( deg2rad * 0.5 * (y1 + y0) );
            double cosTheDx = cosTheta * dx;

            // edge integral
            data[edgeId] = edgeSign * (u[edgeId]*cosTheDx + v[edgeId]*dy);
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

            // ADD CARTESIAN TO DO 
        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__fromVectorFieldEdgeCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]) {

    int numFailures = 0;
    int ier;

    double coef = 0;
    if ((*self)->grid->degrees) {
        coef = 1.0;
    }
    double deg2rad = M_PI / 180.;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            ier = mnt_grid_getPoints(&(*self)->grid, icell, ie, point0, point1);
            numFailures += ier;

            double x0 = point0[LON_INDEX] * deg2rad;
            double x1 = point1[LON_INDEX] * deg2rad;
            double dx = x1 - x0;
            double y0 = point0[LAT_INDEX] * deg2rad;
            double y1 = point1[LAT_INDEX] * deg2rad;
            double dy = y1 - y0;

            // length on the surface of the sphere (1 if no spherical)
            double cosTheDx = cos( coef * 0.5 * (y1 + y0) ) * dx;

            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + ie;
            data[k] = u[k]*cosTheDx + v[k]*dy;

            // ADD CARTESIAN COMPUTATION AND SELECT (TO DO)
        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__fromVectorFieldFaceCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]) {

    int numFailures = 0;
    int ier;

    double deg2rad = M_PI / 180.;

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + ie;

            ier = mnt_grid_getPoints(&(*self)->grid, icell, ie, point0, point1);
            numFailures += ier;

            double lam0 = point0[LON_INDEX] * deg2rad;
            double lam1 = point1[LON_INDEX] * deg2rad;
            double lamMid = 0.5*(lam0 + lam1);

            double the0 = point0[LAT_INDEX] * deg2rad;
            double the1 = point1[LAT_INDEX] * deg2rad;
            double theMid = 0.5*(the0 + the1);
           
            double uval = u[k];
            double vval = v[k];

            double a = 1.0;
            double cosLam = cos(lamMid);
            double sinLam = sin(lamMid);
            double cosThe = cos(theMid);
            double sinThe = sin(theMid);

            // Cartesian estimate
            double zMid = a*sinThe;
            double rhoMid = a*cosThe;
            double xMid = rhoMid*cosLam;
            double yMid = rhoMid*sinLam;

            double vx = - uval*sinLam - vval*zMid*cosLam/a;
            double vy = + uval*cosLam - vval*zMid*sinLam/a;
            double vz = vval*rhoMid/a;

            double rho0 = a*cos(the0);
            double rho1 = a*cos(the1);
            double dx = rho1*cos(lam1) - rho0*cos(lam0);
            double dy = rho1*sin(lam1) - rho0*sin(lam0);
            double dz = a*(sin(the1) - sin(the0));

            double dataXYZ = (dy*zMid - dz*yMid)*vx/a
                           + (dz*xMid - dx*zMid)*vy/a
                           + (dx*yMid - dy*xMid)*vz/a;

            // estimate in spherical coordinates
            double cosTheDLam = cosThe*(lam1 - lam0);
            double dThe = the1 - the0;
            double dataSph = u[k]*dThe - v[k]*cosTheDLam;

            // choose the spherical estimate except when one of the points is at the pole
            const double eps = 10 * std::numeric_limits<double>::epsilon();
            if ((*self)->grid->degrees && 
                (std::fabs(std::fabs(point0[LAT_INDEX]) - 90.) < eps || std::fabs(std::fabs(point1[LAT_INDEX]) - 90.) < eps)) {
                // use estimate obtained from Cartesian coordinates when one of the points is at the poles
                data[k] = dataXYZ;
            }
            else {
                // use estimate from spherical coordinates
                data[k] = dataSph;
            }

        }
    }

    return numFailures;
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

            // TO CHECK: do we want dx, dy in rad or degrees? We will need to be consistent between the
            // toVector and fromVector calls

            // length on the surface of the sphere (1 if no spherical)
            double cosTheDx = cos( deg2rad * 0.5 * (y1 + y0) ) * dx;

            // our convention is to have the extensive fluxes pointing in the positive, logical direction.
            // This requires flipping the sign for iedge = 0 and 2.
            double sign = -2*((iedge + 1) % 2) + 1;

            double len_sq = cosTheDx * cosTheDx + dy * dy;
            if (len_sq <= 0) {
                // happens at the pole
                // edgeData and faceData should be zero there and u, v are ill defined there
                // set to zero (?)
                u[edgeId] = 0;
                v[edgeId] = 0;
            }
            else {
                // normal case
                u[edgeId] = edgeSign * (cosTheDx * edgeData[edgeId] + sign * dy * faceData[edgeId] ) / len_sq;
                v[edgeId] = edgeSign * (dy * edgeData[edgeId] - sign * cosTheDx * faceData[edgeId] ) / len_sq;
            }

        }
    }

    return numFailures;
}

int mnt_extensivefieldadaptor__toVectorFieldCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                           const double* edgeData,
                                                           const double* faceData,
                                                           double* u, double* v) {
    std::string msg;
    int ier, numFailures = 0;

    double deg2rad = 0;
    if ((*self)->grid->degrees) {
        deg2rad = M_PI/180.0;
    }

    double point0[3];
    double point1[3];
    for (vtkIdType icell = 0; icell < (vtkIdType) (*self)->numCells; ++icell) {

        for (int iedge = 0; iedge < MNT_NUM_EDGES_PER_QUAD; ++iedge) {

            ier = mnt_grid_getPoints(&(*self)->grid, icell, iedge, point0, point1);
            if (ier != 0) numFailures++;

            double x0 = point0[LON_INDEX];
            double x1 = point1[LON_INDEX];
            double dx = x1 - x0;
            double y0 = point0[LAT_INDEX];
            double y1 = point1[LAT_INDEX];
            double dy = y1 - y0;

            // TO CHECK: do we want dx, dy in rad or degrees? We will need to be consistent between the
            // toVector and fromVector calls

            // length on the surface of the sphere (1 if no spherical)
            double cosTheDx = cos( deg2rad * 0.5 * (y1 + y0) ) * dx;

            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + iedge;

            double len_sq = cosTheDx * cosTheDx + dy * dy;
            if (len_sq <= 0) {
                // happens at the pole
                // edgeData and faceData should be zero there and u, v are ill defined there
                // set to zero (?)
                u[k] = 0;
                v[k] = 0;
            }
            else {
                // TO CHECK: signs???
                // normal case
                u[k] = (cosTheDx * edgeData[k] + dy * faceData[k] ) / len_sq;
                v[k] = (dy * edgeData[k] - cosTheDx * faceData[k] ) / len_sq;
            }

        }
    }

    return numFailures;
}
