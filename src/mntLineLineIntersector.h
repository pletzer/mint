#ifndef MNT_LINE_LINE_INTERSECTOR
#define MNT_LINE_LINE_INTERSECTOR

#include "mntMatMxN.h"
#include "mntVecN.h"
#include <iostream>
#include <algorithm> // swap()

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

        // start with filling the vectors for the 2D case
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

        double dx = px - qx;
        double dy = py - qy;
        double dz = pz - qz;

        this->solTimesDet[0] = (dz*vy - dy*vz)*wx + (dx*vz + dz*vx)*wy + (dy*vx - dx*vy)*wz;
        this->solTimesDet[1] = (dz*uy - dy*uz)*wx + (dx*uz + dz*ux)*wy + (dy*ux - dx*uy)*wz;
        // we don't need to know the value of the solution in the direction perpendicular 
        // to u and v, so no this->solTimesDet[2]
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

        const Vec3 u = this->p1 - this->p0;
        double uSquare = dot(u, u);
        // lambda @ q0
        double lm0 = dot(this->q0 - this->p0, u)/uSquare;
        // lambda @ q1
        double lm1 = dot(this->q1 - this->p0, u)/uSquare;
        // make sure lambda is in [0, 1]
        this->lamBeg = std::min(std::max(lm0, 0.0), 1.0);
        this->lamEnd = std::min(std::max(lm1, 0.0), 1.0);
        this->lamBegRaw = lm0;
        if (this->lamEnd < this->lamBeg) {
            std::swap(this->lamEnd, this->lamBeg);
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
            // the two lines are on top of each other

            // determinant is zero, p1 - p0 and q1 - q0 are parallel
            this->computeBegEndParamCoords();

            if (std::abs(this->lamEnd - this->lamBeg) > tol) {
                // the two lines overlap

                // compute the distance
                Vec3 p2q = this->q0 - ((1.0 - this->lamBegRaw)*this->p0 + this->lamBegRaw*this->p1);
                double dist = sqrt(dot(p2q, p2q));
                if (dist < tol) {
                    // the lines are touching
                    return true;
                }
            }

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
    double lamBegRaw;
    double det;

};

#endif // MNT_LINE_LINE_INTERSECTOR
