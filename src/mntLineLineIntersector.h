#ifndef MNT_LINE_LINE_INTERSECTOR
#define MNT_LINE_LINE_INTERSECTOR

#include "mntMatMxN.h"
#include "mntVecN.h"
#include <iostream>

struct LineLineIntersector {

    /**
     * Constructor
     */
    LineLineIntersector() {
    }

    /**
     * Set points
     * @param ndim either 2 (2D) or 3 (3D), 
     * @param p0 starting point of first line
     * @param p1 end point of first line
     * @param q0 starting point of second line
     * @param q1 end point of second line
     */
    void setPoints(int ndim, const double p0[], const double p1[], 
                   const double q0[], const double q1[]) {

        // start with filling for the 2D case
        double ux = p1[0] - p0[0];
        double uy = p1[1] - p0[1];
        double uz = 0.;

        double vx = q1[0] - q0[0];
        double vy = q1[1] - q0[1];
        double vz = 0.;

        double wx = 0.;
        double wy = 0.;
        double wz = 1.;

        this->p0[0] = p0[0];
        this->p0[1] = p0[1];
        this->p0[2] = 0.;

        this->p1[0] = p1[0];
        this->p1[1] = p1[1];
        this->p1[2] = 0.;

        this->q0[0] = q0[0];
        this->q0[1] = q0[1];
        this->q0[2] = 0.;

        this->q1[0] = q1[0];
        this->q1[1] = q1[1];
        this->q1[2] = 0.;

        if (ndim == 3) {

            uz = p1[2] - p0[2];
            vz = q1[2] - q0[2];

            // w is perpendicular to both u and v (w = u x v)
            wx = uy*vz - uz*vy;
            wy = uz*vx - ux*vz;
            wz = ux*vy - uy*vx;

            this->p0[2] = p0[2];
            this->p1[2] = p1[2];
            this->q0[2] = q0[2];
            this->q1[2] = q1[2];

        }

        this->det = (uy*vz - uz*vy)*wx + (uz*vx - ux*vz)*wy + (ux*vy - uy*vx)*wz;

        double px = this->p0[0];
        double py = this->p0[1];
        double pz = this->p0[2];

        double qx = this->q0[0];
        double qy = this->q0[1];
        double qz = this->q0[2];

        this->solTimesDet[0] = ((pz - qz)*vy + (qy - py)*vz)*wx + ((qz - pz)*vx + (px - qx)*vz)*wy + ((py - qy)*vx + (qx - px)*vy)*wz;
        this->solTimesDet[1] = ((pz - qz)*uy + (qy - py)*uz)*wx + ((qz - pz)*ux + (px - qx)*uz)*wy + ((py - qy)*ux + (qx - px)*uy)*wz;
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

        const Vec3 dp = this->p1 - this->p0;

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

    Vec2 solTimesDet;

    Vec3 p0;
    Vec3 p1;
    Vec3 q0;
    Vec3 q1;

    double lamBeg;
    double lamEnd;
    double det;

};

#endif // MNT_LINE_LINE_INTERSECTOR
