#include <mntMatMxN.h>
#include <mntVecN.h>
#include <mntLineLineIntersector.h>
#include <iostream>
#include <limits>

#ifndef MNT_LINE_TRIANGLE_INTERSECTOR
#define MNT_LINE_TRIANGLE_INTERSECTOR

#define BAD std::numeric_limits<double>::max()

/**
 * Get the position of p in the triangle parameter space of (q0, q1, q2)
 *
 * @param q0 first triangle vertex (3-element vector)
 * @param q1 second triangle vertex (3-element vector)
 * @param q2 third triangle vertex (3-element vector)
 * @param p target point (3-element vector)
 * @return (xi1, xi2) vector
 * 
 * @note assumes p to be in the plane of (q0, q1, q2)
 */
Vec2
getTriangleParamLocation(const Vec3& q0, 
                         const Vec3& q1, 
                         const Vec3& q2, 
                         const Vec3& p);


struct LineTriangleIntersector {

    /**
     * Constructor
     */
    LineTriangleIntersector();

    /**
     * Set points 
     * @param p0 starting point of line
     * @param p1 end point of line
     * @param q0 first triangle point
     * @param q1 second triangle point
     * @param q2 third triangle point
     */
    void setPoints(const double p0[], 
                   const double p1[], 
                   const double q0[], 
                   const double q1[],
                   const double q2[]);


    /**
     * Get the determinant
     * @return value
     */
    double getDet() const;


    /**
     * Chck if the line runs tangential to the triangle
     * @para tol small tolerance 
     * @return value
     */
    bool isSingular(double tol) const;


    /**
     * Compute the begin/end parametric coordinates
    */
    void computeBegEndParamCoords();

    /**
     * Get the begin/end points of overlap
     * @return pair of points
     */
    const std::pair< Vec3, Vec3 > getBegEndPoints() const;

    /**
     * Get the begin/end parametric coordinates of overlap
     * @return pair
     */
    const std::pair< double, double > getBegEndParamCoords() const;


    /**
     * Check if there is a solution
     * @param tol tolerance
     * @return True if there is one or more solutions
     */
    bool hasSolution(double tol);


    /**
     * Get the solution
     * @return solution
     */
    const Vec3 getSolution();


    Mat3x3 mat;
    Mat3x3 invMatTimesDet;

    Vec3 invMatTimesDetDotRhs;
    Vec3 rhs;
    Vec3 solTimesDet;
    Vec3 p0;
    Vec3 p1;
    Vec3 q0;
    Vec3 q1;
    Vec3 q2;

    double lamBeg;
    double lamEnd;
    double det;

};

#endif // MNT_LINE_TRIANGLE_INTERSECTOR
