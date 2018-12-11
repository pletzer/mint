#include <MvMatrix.h>
#include <MvVector.h>

#ifndef MNT_LINE_TRIANGLE_INTERSECTOR
#define MNT_LINE_TRIANGLE_INTERSECTOR

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
        const Vector<double> dp = this->p1 - this->p0;
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
