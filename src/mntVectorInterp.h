#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVecN.h>

#ifndef MNT_VECTOR_INTERP
#define MNT_VECTOR_INTERP

// index of each vertex
#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

struct VectorInterp_t {
    vmtCellLocator* locator;
    Grid_t* grid;
    Vec3 targetPoint;
    Vec3 pcoords;
    vtkIdType cellId;
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
int mnt_vectorinterp_set(VectorInterp_t** self, Grid_t* grid, vmtCellLocator* locator);

/**
 * Find the target point in the grid
 * @param self instance of VectorInterp_t
 * @param targetPoint target point
 * @param tol2 tolerance
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_find(VectorInterp_t** self, const double targerPoint[], double tol2);

/**
 * Get the contravarian vector components at the parametric cell location
 * @param self instance of VectorInterp_t
 * @param data edge centred data
 * @param vector components (output)
 * @return error code (0 = OK)
 */
extern "C"
int mnt_vectorinterp_getVector(VectorInterp_t** self, const double data[], double vector[]);


#endif // MNT_VECTOR_INTERP
