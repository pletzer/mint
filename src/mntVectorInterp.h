#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVecN.h>

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
extern "C" 
int mnt_vectorinterp_new(VectorInterp_t** self);

/**
 * Destructor
 * @param self instance of VectorInterp_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_del(VectorInterp_t** self);

/**
 * Set the grid and cell locator
 * @param self instance of VectorInterp_t
 * @param grid grid (borrowed reference)
 * @param locator locator (borrowed reference)
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_setGrid(VectorInterp_t** self, Grid_t* grid);

/**
 * Set the grid and cell locator
 * @param self instance of VectorInterp_t
 * @param locator locator (borrowed reference)
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_setLocator(VectorInterp_t** self, vmtCellLocator* locator);

/**
 * Build cell locator
 * @param self instance of VectorInterp_t
 * @param locator locator (borrowed reference)
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_buildLocator(VectorInterp_t** self, int numCellsPerBucket, double periodX);

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
extern "C"
int mnt_vectorinterp_findPoints(VectorInterp_t** self, std::size_t numPoints, 
                                const double targetPoints[], double tol2);

/**
 * Get the edge vectors at given target points
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size numCells*4
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
extern "C"
int mnt_vectorinterp_getEdgeVectors(VectorInterp_t** self,
                                    const double data[], double vectors[]);

/**
 * Get the face vectors at given target points
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size numCells*4
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
extern "C"
int mnt_vectorinterp_getFaceVectors(VectorInterp_t** self,
                                    const double data[], double vectors[]);

/* private */

inline double crossDotZHat(const Vec3& a, const Vec3& b) {
    return a[0]*b[1] - a[1]*b[0];
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
            std::cerr << "Warning: bad cell " << cellId << " vertices: " 
                      << v0 << ',' << v1 << ',' << v2 << ',' << v3 << '\n';
            ier = 1;
        }

        // cotangent vectors obtained by finite differencing and linearly interpolating
        // in the other direction
        drdXsi = ate*a + eta*c;
        drdEta = isx*d + xsi*b;

        return ier;
}
#endif // MNT_VECTOR_INTERP
