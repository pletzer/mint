#include <MvMatrix.h>
#include <MvVector.h>
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
Vector<double>
getTriangleParamLocation(const Vector<double>& q0, 
                         const Vector<double>& q1, 
                         const Vector<double>& q2, 
                         const Vector<double>& p);

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
    const std::pair< Vector<double>, Vector<double> > getBegEndPoints() const;

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
    const Vector<double> getSolution();


    ColMat<double> mat;
    ColMat<double> invMatTimesDet;

    Vector<double> invMatTimesDetDotRhs;
    Vector<double> rhs;
    Vector<double> solTimesDet;
    Vector<double> p0;
    Vector<double> p1;
    Vector<double> q0;
    Vector<double> q1;
    Vector<double> q2;

    double lamBeg;
    double lamEnd;
    double det;

};

#endif // MNT_LINE_TRIANGLE_INTERSECTOR
