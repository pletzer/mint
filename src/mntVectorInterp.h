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
 * Get the vectors at the target points
 * @param self instance of VectorInterp_t
 * @param data edge data, array of size numCells*4
 * @param vectors array of output vectors, size numPoints*3
 * @return error code (0 = OK)
 * @note call this after mnt_vectorinterp_findPoints. The returned vectors
 *       will not be touched if the point falls out of the domain.
 */
extern "C"
int mnt_vectorinterp_getVectors(VectorInterp_t** self,
                                const double data[], double vectors[]);

#endif // MNT_VECTOR_INTERP
