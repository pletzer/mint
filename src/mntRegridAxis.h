#include <vector>
#include <map>
#include <string>
#include <vtkCardinalSpline.h>

#ifndef MNT_REGRID_AXIS
#define MNT_REGRID_AXIS

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct RegridAxis_t {

    // spline object
    vtkCardinalSpline* spline;

    // left/right values in index space
    double tLo, tHi;
    
    // used for cell interpolation
    double tA, tB;
    int indicesA[2], indicesB[2];

    int numCellWeights;
    int numValues;

    // whether of note the axis is increasing monotonically
    bool increasing;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridaxis_new(RegridAxis_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridaxis_del(RegridAxis_t** self);


/**
 * Build the regridder
 * @param numValues number of axis values
 * @param srcValues axis values
 * @return error code (0 is OK, any other value indicates that srcValues are not monotonically increasing/decreasing)
 */
extern "C"
int mnt_regridaxis_build(RegridAxis_t** self, int numValues, const double srcValues[]);

/**
 * Get the point interpolation weights
 * @param target target axis point
 * @param indices indices of the weights (output)
 * @param weights returned weights for the above indices (output)
 * @return error code (0 is OK, 1 target is below range, 2 target is above range)
 */
extern "C"
int mnt_regridaxis_getPointWeights(RegridAxis_t** self, double target, int indices[2], double weights[2]);

/**
 * Get the number of cell interpolation weights
 * @param targets low/high target axis points
 * @param numCellWeights number of returned cell interpolation weights (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridaxis_getNumCellWeights(RegridAxis_t** self, const double targets[2], int* numCellWeights);

/**
 * Get the cell interpolation weights
 * @param indices indices of the weights (filled in)
 * @param weights returned weights for the above indices (filled in)
 * @return error code (0 is OK)
 * @note indices and weights must be pre-allocated. Call mnt_regridaxis_getNumCellWeights to determine the number
 *       of weights
 */
extern "C"
int mnt_regridaxis_getCellWeights(RegridAxis_t** self, int indices[], double weights[]);


#endif // MNT_REGRID_AXIS
