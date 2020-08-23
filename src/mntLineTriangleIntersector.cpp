#include <mntLineTriangleIntersector.h>
#include <mntMatMxN.h>
#include <mntVecN.h>
#include <vector>

Vec2
getTriangleParamLocation(const Vec3& q0, 
                         const Vec3& q1, 
                         const Vec3& q2, 
                         const Vec3& p) {

    Mat3x2 mat;
    Vec3 rhs;
    for (size_t i = 0; i < 3; ++i) {
        mat(i, 0) = q1[i] - q0[i];
        mat(i, 1) = q2[i] - q0[i];
        rhs[i] = p[i] - q0[i];
    }
    Mat2x3 matT = transpose(mat);

    Mat2x2 a = dot(matT, mat);

    double det = a(0, 0)*a(1, 1) - a(0, 1)*a(1, 0);
    Mat2x2 aInv;
    aInv(0, 0) = a(1, 1)/det;
    aInv(0, 1) = -a(0, 1)/det;
    aInv(1, 0) = -a(1, 0)/det;
    aInv(1, 1) = a(0, 0)/det;

    return dot(aInv, dot(matT, rhs));
}


LineTriangleIntersector::LineTriangleIntersector() {
    this->lamBeg = BAD;
    this->lamEnd = BAD;
}

void 
LineTriangleIntersector::setPoints(const double p0[],
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


double 
LineTriangleIntersector::getDet() const {
    return this->det;
}


bool 
LineTriangleIntersector::isSingular(double tol) const {
    if (std::abs(this->det) < tol)
        return true;
    return false;
}


void 
LineTriangleIntersector::computeBegEndParamCoords() {

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
    Vec2 pxi0 = getTriangleParamLocation(this->q0, this->q1, this->q2, this->p0);
    Vec2 pxi1 = getTriangleParamLocation(this->q0, this->q1, this->q2, this->p1);

    const double tol = 1.e-10;

    // collect intersection points
    std::vector<double> lamVals;

    if (pxi0[0] > 0. - tol && pxi0[1] < 1. - pxi0[0] + tol) {
        // add this point if it is in the triangle
        lamVals.push_back(0.0);
    }
    if (pxi1[0] > 0. - tol && pxi1[1] < 1. - pxi1[0] + tol) {
        // add this point if it is in the triangle
        lamVals.push_back(1.0);
    }

    llA.setPoints(2, &pxi0[0], &pxi1[0], qxi0, qxi1);
    llB.setPoints(2, &pxi0[0], &pxi1[0], qxi1, qxi2);
    llC.setPoints(2, &pxi0[0], &pxi1[0], qxi2, qxi0);

    Vec2 sol;

    if (llA.hasSolution(tol)) {
        sol = llA.getSolution();
        if (sol[0] > 0. - tol && sol[0] < 1.0 + tol && 
            sol[1] > 0. - tol && sol[1] < 1.0 + tol) {
            lamVals.push_back(sol[0]);
        }
    }
    if (llB.hasSolution(tol)) {
        sol = llB.getSolution();
        if (sol[0] > 0. - tol && sol[0] < 1.0 + tol && 
            sol[1] > 0. - tol && sol[1] < 1.0 + tol) {
            lamVals.push_back(sol[0]);
        }
    }
    if (llC.hasSolution(tol)) {
        sol = llC.getSolution();
        if (sol[0] > 0. - tol && sol[0] < 1.0 + tol && 
            sol[1] > 0. - tol && sol[1] < 1.0 + tol) {
            lamVals.push_back(sol[0]);
        }
    }

    if (lamVals.size() >= 2) {
        // sort
        std::sort(lamVals.begin(), lamVals.end());

        // set
        this->lamBeg = lamVals[0];
        this->lamEnd = lamVals[lamVals.size() - 1];
    }
}


const std::pair< Vec3, Vec3 > 
LineTriangleIntersector::getBegEndPoints() const {
    Vec3 dp = this->p1 - this->p0;
    std::pair< Vec3, Vec3 > p(this->p0 + this->lamBeg*dp,
                              this->p0 + this->lamEnd*dp);
    return p;
}


const std::pair< double, double > 
LineTriangleIntersector::getBegEndParamCoords() const {
       return std::pair< double, double >(this->lamBeg, this->lamEnd);
    }


bool 
LineTriangleIntersector::hasSolution(double tol) {

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


const Vec3
LineTriangleIntersector::getSolution() {
    Vec3 res = this->solTimesDet;
    res /= this->det;
    return res;
}


