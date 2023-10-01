#include "mntLIBRARY_API.h"
#include <mntGlobal.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVecN.h>
#include <sstream> // std::stringstream
#include "mntLogger.h"

#ifndef MNT_VECTOR_INTERP
#define MNT_VECTOR_INTERP


struct VectorInterp_t {

    // grid
    Grid_t* grid;

    // found cell Ids, one for each target point
    std::vector<vtkIdType> cellIds;

    // parametric coordinates for the quad, one per
    // target point
    std::vector<Vec3> pcoords;

    // cell locator, either a borrowed pointer or
    // or owned by the class
    vmtCellLocator* locator;

    bool ownsLocator;
    bool calledFindPoints;
};

/**
 * Constructor
 * @param self instance of VectorInterp_t
 * @return error code (0 = OK)
 */
LIBRARY_API 
int mnt_vectorinterp_new(VectorInterp_t** self);

/**
 * Destructor
 * @param self instance of VectorInterp_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_del(VectorInterp_t** self);

/**
 * Set the grid
 * @param self instance of VectorInterp_t
 * @param grid grid (borrowed reference)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_setGrid(VectorInterp_t** self, Grid_t* grid);

/**
 * Set the grid cell locator
 * @param self instance of VectorInterp_t
 * @param locator locator (borrowed reference)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_setLocator(VectorInterp_t** self, vmtCellLocator* locator);

/**
 * Build the grid cell locator
 * @param self instance of VectorInterp_t
 * @param numCellsPerBucket number of cells per bucket. The smaller the faster the cell search. However, 
 *                          small values may casue problems, we recommend about 100 or more
 * @param periodX period length, use 0 if non-periodic in the first coordinate
 * @param enableFolding whether (1) or not (0) |latitude| > 90 deg values should be folded back into the domain
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_buildLocator(VectorInterp_t** self, int numCellsPerBucket, double periodX, int enableFolding);

/**
 * Find target points
 * @param self instance of VectorInterp_t
 * @param numPoints number of points
 * @param targetPoints array of target points, size numPoints*3
 * @param tol2 tolerance
 * @return error code (0 = OK)
 * @note the error code is the number of points for which the parametric coordinates could 
 *       not be found
 */
LIBRARY_API
int mnt_vectorinterp_findPoints(VectorInterp_t** self, std::size_t numPoints, 
                                const double targetPoints[], double tol2);

/**
 * Get the edge vectors at given target points from cell by cell data
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size numCells*MNT_NUM_EDGES_PER_QUAD. The units must be consistent
 *             with the coordinates.
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp__getEdgeVectorsFromCellByCellData(VectorInterp_t** self,
                                                      const double data[],
                                                      double vectors[]);

/**
 * Get the face vectors at given target points from cell by cell data
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size numCells*MNT_NUM_EDGES_PER_QUAD. The units must be consistent
 *             with the coordinates.
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp__getFaceVectorsFromCellByCellData(VectorInterp_t** self,
                                                      const double data[],
                                                      double vectors[]);
/**
 * Get the edge vectors at given target points from unique edge Id data
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size number of edges. The units must be
 *             consistent with the coordinates.
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp__getEdgeVectorsFromUniqueEdgeData(VectorInterp_t** self,
                                                      const double data[],
                                                      double vectors[]);

/**
 * Get the face vectors at given target points from unique edge Id data
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size number of edges. The units must be
 *             consistent with the coordinates.
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp__getFaceVectorsFromUniqueEdgeData(VectorInterp_t** self,
                                                      const double data[],
                                                      double vectors[]);

/**
 * Get the edge vectors at given target points
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data
 * @param placement MNT_CELL_BY_CELL_DATA if data are cell by cell 
 *                  (size of array is numCells * MNT_NUM_EDGES_PER_QUAD),
 *                  otherwise data are assumed to be on unique edges 
 *                  (size is numEdges)
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp_getEdgeVectors(VectorInterp_t** self,
                                    const double data[], int placement,
                                    double vectors[]);

/**
 * Get the Face vectors at given target points
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data
 * @param placement MNT_CELL_BY_CELL_DATA if data are cell by cell
 *                  (size numCells * MNT_NUM_EDGES_PER_QUAD),
 *                  otherwise data are assumed to be on unique edges 
 *                  (size is numEdges)
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
LIBRARY_API
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self,
                                    const double data[], int placement,
                                    double vectors[]);

/**
 * Get the vectors at the mid edge locations from the extensive field in Cartesian coordinates
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data, size depends on placement (see below)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. 
 *                  If placement == MNT_CELL_BY_CELL_DATA then data size is 
 *                  num cells * num edges per cell, if MNT_UNIQUE_EDGE_DATA then
 *                  data size should be num edges.
 * @param u x-component of the output vectors, size numEdges
 * @param v y-component of the output vectors, size numEdges
 * @param fs function space, either MNT_FUNC_SPACE_W1 or MNT_FUNC_SPACE_W2
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_getVectorsOnEdges(VectorInterp_t** self,
                                        const double data[],
                                        int placement,
                                        double u[], double v[],
                                        int fs);


/**
 * Get the vectors at the mid edge locations from the extensive field in spherical coordinates
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data, size depends on placement (see below)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. 
 *                  If placement == MNT_CELL_BY_CELL_DATA then data size is 
 *                  num cells * num edges per cell, if MNT_UNIQUE_EDGE_DATA then
 *                  data size should be num edges.
 * @param u x-component of the output vectors, size numEdges
 * @param v y-component of the output vectors, size numEdges
 * @param fs function space, either MNT_FUNC_SPACE_W1 or MNT_FUNC_SPACE_W2
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_getVectorsOnEdgesSpherical(VectorInterp_t** self,
                                            const double data[],
                                            int placement,
                                            double u[], double v[],
                                            int fs);

/* private */

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    Vec3 res;
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}

inline double crossDotZHat(const Vec3& a, const Vec3& b) {
    return a[0]*b[1] - a[1]*b[0];
}

inline double crossDot(const Vec3& a, const Vec3& b, const Vec3& c) {
    return c[0]*(a[1]*b[2] - a[2]*b[1]) + 
           c[1]*(a[2]*b[0] - a[0]*b[2]) + 
           c[2]*(a[0]*b[1] - a[1]*b[0]);
}

inline Vec3 cartesianFromRadians(const Vec3& p) {
    Vec3 res;
    double lam = p[LON_INDEX];
    double the = p[LAT_INDEX];
    res[0] = cos(the) * cos(lam);
    res[1] = cos(the) * sin(lam);
    res[2] = sin(the);
    return res;
}


/*
 * Get the edge vector at some cell target location
 * @param v0 1st vertex coordinates of the quad
 * @param v1 2nd vertex coordinates of the quad
 * @param v2 3rd vertex coordinates of the quad
 * @param v3 4th vertex coordinates of the quad
 * @param xsi 1st parametric cooridnate of the target
 * @param eta 2nd parametric cooridnate of the target
 * @param edgeValues edge integrated values, assumes the edges to point in the positive xsi, eta directions
 * @param vector output vector at the xsi, eta location within the v0, v1, v2, v3 cell
 * @return 0 if no error
 */
inline int mnt_vectorinterp__getEdgeVector(const Vec3& v0, 
                                           const Vec3& v1,
                                           const Vec3& v2,
                                           const Vec3& v3,
                                           double xsi,
                                           double eta,
                                           const double edgeValues[], 
                                           double vector[]) {

    double isx = 1.0 - xsi;
    double ate = 1.0 - eta;

    Vec3 a = v1 - v0;
    Vec3 b = v2 - v1;
    Vec3 c = v2 - v3;
    Vec3 d = v3 - v0;

    // cotangent vectors obtained by finite differencing and linearly interpolating
    // in the other direction
    Vec3 drdXsi = ate*a + eta*c;
    Vec3 drdEta = isx*d + xsi*b;

    Vec3 area = cross(drdXsi, drdEta);
    Vec3 normal = area / sqrt(dot(area, area));
    double jac = dot(area, normal);
    if (jac <= 0) {
        std::stringstream msg;
        msg << "bad cell: vertices: " <<
                         v0 << ";" << v1 << ";" << v2  << ";" << v3; 
        mntlog::warn(__FILE__, __func__, __LINE__, msg.str());
        return 1;
    }
    Vec3 gradXsi = cross(drdEta, normal) / jac;
    Vec3 gradEta = cross(normal, drdXsi) / jac;

    for (auto i = 0; i < 3; ++i) {
        vector[i] = (ate*edgeValues[0] + eta*edgeValues[2]) * gradXsi[i] + 
                    (isx*edgeValues[3] + xsi*edgeValues[1]) * gradEta[i];
    }

    return 0;
}


/*
 * Get the face vector at some cell target location
 * @param v0 1st vertex coordinates of the quad
 * @param v1 2nd vertex coordinates of the quad
 * @param v2 3rd vertex coordinates of the quad
 * @param v3 4th vertex coordinates of the quad
 * @param xsi 1st parametric cooridnate of the target
 * @param eta 2nd parametric cooridnate of the target
 * @param edgeValues edge integrated values, assumes the edges to point in the positive xsi, eta directions
 * @param vector output vector at the xsi, eta location within the v0, v1, v2, v3 cell
 * @return 0 if no error
 */
inline int mnt_vectorinterp__getFaceVector(const Vec3& v0, 
                                           const Vec3& v1,
                                           const Vec3& v2,
                                           const Vec3& v3,
                                           double xsi,
                                           double eta,
                                           const double edgeValues[], 
                                           double vector[]) {

    double isx = 1.0 - xsi;
    double ate = 1.0 - eta;

    Vec3 a = v1 - v0;
    Vec3 b = v2 - v1;
    Vec3 c = v2 - v3;
    Vec3 d = v3 - v0;

    // cotangent vectors obtained by finite differencing and linearly interpolating
    // in the other direction
    Vec3 drdXsi = ate*a + eta*c;
    Vec3 drdEta = isx*d + xsi*b;

    Vec3 area = cross(drdXsi, drdEta);
    Vec3 normal = area / sqrt(dot(area, area));
    double jac = dot(area, normal);
    if (jac <= 0) {
        std::stringstream msg;
        msg << "bad cell: vertices: " <<
                         v0 << ";" << v1 << ";" << v2  << ";" << v3; 
        mntlog::warn(__FILE__, __func__, __LINE__, msg.str());
        return 1;
    }
    Vec3 gradXsi = cross(drdEta, normal) / jac;
    Vec3 gradEta = cross(normal, drdXsi) / jac;
    Vec3 nXGradXsi = cross(normal, gradXsi);
    Vec3 gradEtaXn = cross(gradEta, normal);

    for (auto i = 0; i < 3; ++i) {
        vector[i] = (ate*edgeValues[0] + eta*edgeValues[2]) * nXGradXsi[i] + 
                    (isx*edgeValues[3] + xsi*edgeValues[1]) * gradEtaXn[i];
    }

    return 0;
}


#endif // MNT_VECTOR_INTERP
