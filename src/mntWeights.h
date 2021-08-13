#include <mntVecN.h>

const double QUAD_COUNTERCLOCKWISE_EDGE_DIRECTION[] = {
    0., 0., 0.,   1., 0., 0., // south
    1., 0., 0.,   1., 1., 0., // east
    1., 1., 0.,   0., 1., 0., // north
    0., 1., 0.,   0., 0., 0., // west
};

const double QUAD_POSITIVEXI_EDGE_DIRECTION[] = {
    0., 0., 0.,   1., 0., 0., // south
    1., 0., 0.,   1., 1., 0., // east
    0., 1., 0.,   1., 1., 0., // north
    0., 0., 0.,   0., 1., 0.  // west
};

/**
 * Compute the interpolation weight between a source cell edge and a destination line segment
 * @param srcXi0 starting point of src edge
 * @param srcXi1 end point of the src edge
 * @param xia starting point of target line
 * @param xib end point of target line
 * @return interpolation weight
 */
double computeWeight(const double srcXi0[], const double srcXi1[], const Vec3& xia, const Vec3& xib);
