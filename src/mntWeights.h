#include <mntVecN.h>

/**
 * Compute the interpolation weight between a source cell edge and a destination line segment
 * @param srcXi0 starting point of src edge
 * @param srcXi1 end point of the src edge
 * @param dstXi0 start point of target line
 * @param dstXi1 end point of target line
 * @return interpolation weight
 */
double computeWeight(const double srcXi0[], const double srcXi1[], const Vec3& xia, const Vec3& xib);
