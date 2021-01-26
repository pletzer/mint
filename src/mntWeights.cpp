#include <mntVecN.h>

/**
 * Compute the interpolation weight between a source cell edge and a destination line segment
 * @param srcXi0 starting point of src edge
 * @param srcXi1 end point of the src edge
 * @param dstXi0 start point of target line
 * @param dstXi1 end point of target line
 * @return interpolation weight
 */
double computeWeight(const double srcXi0[], const double srcXi1[],
                     const Vec3& xia, const Vec3& xib) {

    double weight = 1.0;
    double sgn = 0.0;

    // iterate over the two dimensions of the quad
    for (size_t d = 0; d < 2; ++d) { // 2d 

        double xiM = 0.5*(xia[d] + xib[d]);
        double dxi = xib[d] - xia[d];

        // mid point of edge in parameter space
        double xm = 0.5*(srcXi0[d] + srcXi1[d]);
        sgn += srcXi1[d] - srcXi0[d]; // should be either -1, 0 or 1

        // use Lagrange interpolation to evaluate the basis function integral for
        // any of the 3 possible x values in {0, 0.5, 1}. This formula will make 
        // it easier to extend the code to 3d
        double xm00 = xm;
        double xm05 = xm - 0.5;
        double xm10 = xm - 1.0;
        double lag00 = + 2.0 * xm05 * xm10;
        double lag05 = - 4.0 * xm00 * xm10;
        double lag10 = + 2.0 * xm00 * xm05;

        weight *= (1.0 - xiM)*lag00 + dxi*lag05 + xiM*lag10;
    }

    return sgn * weight;
}
