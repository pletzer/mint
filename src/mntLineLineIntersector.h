#ifndef MNT_LINE_LINE_INTERSECTOR
#define MNT_LINE_LINE_INTERSECTOR

#include "mntMatMxN.h"
#include "mntVecN.h"

struct LineLineIntersector {

    /**
     * Constructor
     */
    LineLineIntersector() {
    }

    /**
     * Set points 
     * @param p0 starting point of first line
     * @param p1 end point of first line
     * @param q0 starting point of second line
     * @param q1 end point of second line
     */
    void setPoints(const double p0[], 
                   const double p1[], 
                   const double q0[], 
                   const double q1[]) {
        for (size_t i = 0; i < 2; ++i) {
            this->p0[i] = p0[i];
            this->p1[i] = p1[i];
            this->q0[i] = q0[i];
            this->q1[i] = q1[i];
            this->rhs[i] = q0[i] - p0[i];
            this->mat(i, 0) = p1[i] - p0[i];
            this->mat(i, 1) = q0[i] - q1[i];
        }
        this->invMatTimesDet(0, 0) = this->mat(1, 1);
        this->invMatTimesDet(1, 1) = this->mat(0, 0);
        this->invMatTimesDet(0, 1) = -this->mat(0, 1);
        this->invMatTimesDet(1, 0) = -this->mat(1, 0);
        this->solTimesDet = dot(this->invMatTimesDet, this->rhs);
        this->det = this->mat(0, 0)*this->mat(1, 1) - this->mat(0, 1)*this->mat(1, 0);
    }

    /**
     * Get the determinant
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
        const Vec2 dp = this->p1 - this->p0;
        double dp2 = dot(dp, dp);
        // lambda @ q0
        double lm0 = dot(this->q0 - this->p0, dp)/dp2;
        // lambda @ q1
        double lm1 = dot(this->q1 - this->p0, dp)/dp2;

        this->lamBeg = std::min(std::max(lm0, 0.0), 1.0);
        this->lamEnd = std::max(std::min(lm1, 1.0), 0.0);
        if (this->lamEnd < this->lamBeg) {
            // switch
            double lbeg = this->lamEnd;
            double lend = this->lamBeg;
            this->lamBeg = lbeg;
            this->lamEnd = lend;
        }
    }


    /**
     * Get the begin/end points of overlap
     * @return pair of points
     */
    const std::pair< Vec2, Vec2 > getBegEndPoints() const {
        Vec2 dp = this->p1 - this->p0;
        std::pair< Vec2, Vec2 > p(this->p0 + this->lamBeg*dp, 
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
            // determinant is zero, p1 - p0 and q1 - q0 are on
            // the same ray
            this->computeBegEndParamCoords();

            if (std::abs(this->lamEnd - this->lamBeg) > tol)
                return true;
        }

        return false;
    }


    /**
     * Get the solution
     * @return solution
     */
    const Vec2 getSolution() {
        Vec2 res = this->solTimesDet;
        res /= this->det;
        return res; 
    }

    Mat2x2 mat;
    Mat2x2 invMatTimesDet;

    Vec2 invMatTimesDetDotRhs;
    Vec2 rhs;
    Vec2 solTimesDet;
    Vec2 p0;
    Vec2 p1;
    Vec2 q0;
    Vec2 q1;

    double lamBeg;
    double lamEnd;
    double det;

};

#endif // MNT_LINE_LINE_INTERSECTOR
