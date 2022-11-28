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
 * @param locator locator (borrowed reference)
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
 * Get the edge vectors at the mid edge locations from the extensive field
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data, size depends on placement (see below)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. 
 *                  If placement == MNT_CELL_BY_CELL_DATA then data size is 
 *                  num cells * num edges per cell, if MNT_UNIQUE_EDGE_DATA then
 *                  data size should be num edges.
 * @param u x-component of the output vectors, size numEdges
 * @param v y-component of the output vectors, size numEdges
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_getEdgeVectorsOnEdges(VectorInterp_t** self,
                                            const double data[],
                                            int placement,
                                            double u[], double v[]);

/**
 * Get the face vectors at the mid edge locations from the extensive field
 * @param self instance of VectorInterp_t
 * @param data edge integrated (extensive) data, size depends on placement (see below)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. 
 *                  If placement == MNT_CELL_BY_CELL_DATA then data size is 
 *                  num cells * num edges per cell, if MNT_UNIQUE_EDGE_DATA then
 *                  data size should be num edges.
 * @param u x-component of the output vectors, size numEdges
 * @param v y-component of the output vectors, size numEdges
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vectorinterp_getFaceVectorsOnEdges(VectorInterp_t** self,
                                            const double data[],
                                            int placement,
                                            double u[], double v[]);



/* private */

LIBRARY_API
int mnt_vectorinterp__getEdgeVectorsFromUniqueEdgeDataOnEdges(VectorInterp_t** self,
                                                          const double data[],
                                                          double u[], double v[]);
LIBRARY_API
int mnt_vectorinterp__getEdgeVectorsFromCellByCellDataOnEdges(VectorInterp_t** self,
                                                          const double data[],
                                                          double u[], double v[]);
LIBRARY_API
int mnt_vectorinterp__getFaceVectorsFromUniqueEdgeDataOnEdges(VectorInterp_t** self,
                                                          const double data[],
                                                          double u[], double v[]);
LIBRARY_API
int mnt_vectorinterp__getFaceVectorsFromCellByCellDataOnEdges(VectorInterp_t** self,
                                                          const double data[],
                                                          double u[], double v[]);


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

inline int mnt_vectorinterp__getTangentVectors(VectorInterp_t** self, std::size_t iTargetId,
                                                Vec3& drdXsi, Vec3& drdEta, double& jac) {

        Vec3 v0, v1, v2, v3;
        vtkIdType cellId = (*self)->cellIds[iTargetId];
        int ier = 0;

        // parametric coordinates of the target point 
        double xsi = (*self)->pcoords[iTargetId][0];
        double eta = (*self)->pcoords[iTargetId][1];
        double isx = 1.0 - xsi;
        double ate = 1.0 - eta;

        // get the cell vertices, this should never fail 
        mnt_grid_getPoints(&(*self)->grid, cellId, 0, &v0[0], &v1[0]);
        mnt_grid_getPoints(&(*self)->grid, cellId, 2, &v3[0], &v2[0]);

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
        jac = 0.25*(a013 + a120 + a231 + a302);
        if (jac <= 0) {
            std::stringstream msg;
            msg << "bad cell " << cellId << " vertices: " <<
                            v0 << ";" << v1 << ";" << v2  << ";" << v3; 
            mntlog::warn(__FILE__, __func__, __LINE__, msg.str());
            ier = 1;
        }

        // cotangent vectors obtained by finite differencing and linearly interpolating
        // in the other direction
        drdXsi = ate*a + eta*c;
        drdEta = isx*d + xsi*b;

        return ier;
}
#endif // MNT_VECTOR_INTERP
