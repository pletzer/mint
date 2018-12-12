#include <MvMatrix.h>
#include <MvVector.h>
#include <mntLineLineIntersector.h>
#include <iostream>

#ifndef MNT_LINE_TRIANGLE_INTERSECTOR
#define MNT_LINE_TRIANGLE_INTERSECTOR

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
                         const Vector<double>& p) {
    // 3x2 matrix
    ColMat<double> mat(3, 2);
    Vector<double> rhs(3);
    for (size_t i = 0; i < 3; ++i) {
        mat(i, 0) = q1[i] - q0[i];
        mat(i, 1) = q2[i] - q0[i];
        rhs[i] = p[i] - q0[i];
    }
    ColMat<double> matT = transpose(mat);
    ColMat<double> a = dot(matT, mat);
    double det = a(0, 0)*a(1, 1) - a(0, 1)*a(1, 0);
    ColMat<double> aInv(2, 2);
    aInv(0, 0) = a(1, 1)/det;
    aInv(0, 1) = -a(0, 1)/det;
    aInv(1, 0) = -a(1, 0)/det;
    aInv(1, 1) = a(0, 0)/det;
    return dot(aInv, matT, rhs);
}

struct LineTriangleIntersector {

    /**
     * Constructor
     */
    LineTriangleIntersector() {
        this->mat.newsize(3, 3);
        this->invMatTimesDet.newsize(3, 3);
        this->invMatTimesDetDotRhs.alloc(3);
        this->rhs.alloc(3);
        this->p0.alloc(3);
        this->p1.alloc(3);
        this->q0.alloc(3);
        this->q1.alloc(3);
        this->q2.alloc(3);
    }

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
                   const double q2[]) {

        for (size_t i = 0; i < 3; ++i) {
            this->p0[i] = p0[i];
            this->p1[i] = p1[i];
            this->q0[i] = q0[i];
            this->q1[i] = q1[i];
            this->q2[i] = q2[i];
            this->rhs[i] = q0[i] - p0[i];
            this->mat(i, 0) = p1[i] - p0[i];
            this->mat(i, 1) = q0[i] - q1[i];
            this->mat(i, 2) = q0[i] - q2[i];
        }

        this->invMatTimesDet(0, 0) = -this->mat(1, 2)*this->mat(2, 1) + this->mat(1, 1)*this->mat(2, 2);
        this->invMatTimesDet(0, 1) = +this->mat(0, 2)*this->mat(2, 1) - this->mat(0, 1)*this->mat(2, 2);
        this->invMatTimesDet(0, 2) = -this->mat(0, 2)*this->mat(1, 1) + this->mat(0, 1)*this->mat(1, 2);
        this->invMatTimesDet(1, 0) = +this->mat(1, 2)*this->mat(2, 0) - this->mat(1, 0)*this->mat(2, 2);
        this->invMatTimesDet(1, 1) = -this->mat(0, 2)*this->mat(2, 0) + this->mat(0, 0)*this->mat(2, 2);
        this->invMatTimesDet(1, 2) = +this->mat(0, 2)*this->mat(1, 0) - this->mat(0, 0)*this->mat(1, 2);
        this->invMatTimesDet(2, 0) = -this->mat(1, 1)*this->mat(2, 0) + this->mat(1, 0)*this->mat(2, 1);
        this->invMatTimesDet(2, 1) = +this->mat(0, 1)*this->mat(2, 0) - this->mat(0, 0)*this->mat(2, 1);
        this->invMatTimesDet(2, 2) = -this->mat(0, 1)*this->mat(1, 0) + this->mat(0, 0)*this->mat(1, 1);

        this->solTimesDet = dot(this->invMatTimesDet, this->rhs);

        this->det = - this->mat(0, 2)*this->mat(1, 1)*this->mat(2, 0)
                    + this->mat(0, 1)*this->mat(1, 2)*this->mat(2, 0) 
                    + this->mat(0, 2)*this->mat(1, 0)*this->mat(2, 1) 
                    - this->mat(0, 0)*this->mat(1, 2)*this->mat(2, 1) 
                    - this->mat(0, 1)*this->mat(1, 0)*this->mat(2, 2) 
                    + this->mat(0, 0)*this->mat(1, 1)*this->mat(2, 2);
    }

    /**
     * Get the determinant of the linear system
     * @return determinant
    */
    double getDet() const {
        return this->det;
    }

    /**
     * Check if there is a solution
     * @param tol tolerance
     * @return True if there is one or more solutions
    */
    bool isSingular(double tol) const {
        if (std::abs(this->det) < tol)
            return true;
        return false;
    }


    /**
     * Compute the begin/end parametric coordinates
    */
    void computeBegEndParamCoords() {

        // intersectors of the line with any of the three edges of the
        // triangle
        LineLineIntersector llA;
        LineLineIntersector llB;
        LineLineIntersector llC;

        // endpoints of the triangle in the parametric space of the triangle
        const double qxi0[] = {0., 0.};
        const double qxi1[] = {1., 0.};
        const double qxi2[] = {0., 1.};

        // endpoints of the line in the parametric space of the triangle
        // are computed using a pseudo-inverse
        Vector<double> pxi0 = getTriangleParamLocation(this->q0, this->q1, this->q2, this->p0);
        Vector<double> pxi1 = getTriangleParamLocation(this->q0, this->q1, this->q2, this->p1);

        llA.setPoints(&pxi0[0], &pxi1[0], qxi0, qxi1);
        llB.setPoints(&pxi0[0], &pxi1[0], qxi1, qxi2);
        llC.setPoints(&pxi0[0], &pxi1[0], qxi2, qxi0);

        // collect intersection points
        std::vector<double> lamVals;
        Vector<double> sol;

        const double tol = 1.e-10;
        if (llA.hasSolution(tol)) {
            sol = llA.getSolution();
            if (sol[0] > 0. - tol && sol[0] < 1.0 + tol) {
                lamVals.push_back(sol[0]);
            }
        }
        if (llB.hasSolution(tol)) {
            sol = llB.getSolution();
            if (sol[0] > 0. - tol && sol[0] < 1.0 + tol) {
                lamVals.push_back(sol[0]);
            }
        }
        if (llC.hasSolution(tol)) {
            sol = llC.getSolution();
            if (sol[0] > 0. - tol && sol[0] < 1.0 + tol) {
                lamVals.push_back(sol[0]);
            }
        }

        // sort
        std::sort(lamVals.begin(), lamVals.end());

        // set
        this->lamBeg = lamVals[0];
        this->lamEnd = lamVals[lamVals.size() - 1];
    }


    /**
     * Get the begin/end points of overlap
     * @return pair of points
     */
    const std::pair< Vector<double>, Vector<double> > getBegEndPoints() const {
        Vector<double> dp = this->p1 - this->p0;
        std::pair< Vector<double>, Vector<double> > p(this->p0 + this->lamBeg*dp, 
                                                      this->p0 + this->lamEnd*dp);
        return p;
    }


    /**
     * Get the begin/end parametric coordinates of overlap
     * @return pair
     */
    const std::pair< double, double > getBegEndParamCoords() const {
       return std::pair< double, double >(this->lamBeg, this->lamEnd);
    }


    /**
     * Check if there is a solution
     * @param tol tolerance
     * @return True if there is one or more solutions
     */
    bool hasSolution(double tol) {

        if (std::abs(this->getDet()) > tol)
            return true;

        if (std::abs(dot(this->solTimesDet, this->solTimesDet)) < tol) {
            // determinant is zero
            // ray and triangle are on the same plane
            this->computeBegEndParamCoords();
            std::cerr << "*** lamBeg = " << lamBeg << ", lamEnd = " << lamEnd << '\n';

            if (std::abs(this->lamEnd - this->lamBeg) > tol)
                return true;
        }

        return false;
    }


    /**
     * Get the solution
     * @return solution
     */
    const Vector<double> getSolution() {
        Vector<double> res = this->solTimesDet;
        res /= this->det;
        return res; 
    }

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
