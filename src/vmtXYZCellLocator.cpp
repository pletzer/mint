#include <vmtXYZCellLocator.h>
#include <mntLogger.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <sstream>
#include <limits>
#include <algorithm>


static inline Vec3 crossProduct(const Vec3& a, const Vec3& b) {
    Vec3 res;
    res[0] = a[1]*b[2] - a[2]*b[1];
    res[1] = a[2]*b[0] - a[0]*b[2];
    res[2] = a[0]*b[1] - a[1]*b[0];
    return res;
}


vmtXYZCellLocator::vmtXYZCellLocator() {
    this->locator = vtkStaticCellLocator::New();
    this->scratchCell = vtkGenericCell::New();
    this->grid = nullptr;
    this->projCentre = Vec3(0.0);
}


vmtXYZCellLocator::~vmtXYZCellLocator() {
    this->locator->Delete();
    this->scratchCell->Delete();
}


void
vmtXYZCellLocator::SetDataSet(vtkUnstructuredGrid* grid) {
    this->grid = grid;
    this->locator->SetDataSet(grid);
}


void
vmtXYZCellLocator::SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) {
    this->locator->SetNumberOfCellsPerNode(avgNumFacesPerBucket);
}


void
vmtXYZCellLocator::BuildLocator() {
    this->locator->BuildLocator();
    this->projCentre = Vec3(0.0);
    this->setProjectionCentre();
}


void
vmtXYZCellLocator::setProjectionCentre() {
    // The cells are projected centrally from projCentre, which must see
    // every cell from the same side: the origin for a sphere (or any
    // surface around it), otherwise a point far off the mesh along its
    // mean normal (e.g. for a flat mesh through the origin)
    vtkIdType numCells = this->grid->GetNumberOfCells();
    if (numCells == 0) {
        return;
    }
    Vec3 verts[4];
    Vec3 meanCentroid(0.0), meanNormal(0.0);
    double minCos = 1.0, maxCos = -1.0;
    double bounds[6];
    this->grid->GetBounds(bounds);
    double diam = std::sqrt((bounds[1] - bounds[0])*(bounds[1] - bounds[0]) +
                            (bounds[3] - bounds[2])*(bounds[3] - bounds[2]) +
                            (bounds[5] - bounds[4])*(bounds[5] - bounds[4]));
    for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
        this->getCellVerts(cellId, verts);
        Vec3 c = 0.25 * (verts[0] + verts[1] + verts[2] + verts[3]);
        Vec3 n = crossProduct(verts[2] - verts[0], verts[3] - verts[1]);
        double nn = std::sqrt(dot(n, n)), cc = std::sqrt(dot(c, c));
        if (nn > 0.) {
            meanNormal += n / nn;
        }
        meanCentroid += c;
        double cosAngle = (nn > 0. && cc > 0.) ? dot(n, c) / (nn * cc) : 0.0;
        minCos = std::min(minCos, cosAngle);
        maxCos = std::max(maxCos, cosAngle);
    }
    // every cell faces the origin, from the same side
    if (minCos > 0.1 || maxCos < -0.1) {
        return;
    }
    meanCentroid = meanCentroid / double(numCells);
    double mn = std::sqrt(dot(meanNormal, meanNormal));
    if (mn < 1.e-10 * numCells) {
        mntlog::warn(__FILE__, __func__, __LINE__,
                     "the cells do not all face the origin and have no mean normal: "
                     "line tracing may be inaccurate");
        return;
    }
    this->projCentre = meanCentroid - (100. * diam / mn) * meanNormal;
}


vtkIdType
vmtXYZCellLocator::FindCell(const double point[3], double tol2, vtkGenericCell *cell,
                             double pcoords[3], double *weights) {
    // vtkStaticCellLocator::FindCell needs a real cell to write into, and
    // its own signature is not const-correct about `point` even though it
    // doesn't modify it -- copy into a local, mutable array. Deliberately
    // using the 5-argument overload (no subId): some VTK versions (e.g.
    // 9.1, as shipped by several Linux distros) only expose that one on
    // vtkStaticCellLocator itself -- declaring any FindCell override hides
    // the base class's other overloads (including the 6-argument one with
    // subId) unless re-exposed with a `using` declaration, which VTK 9.1
    // does not do here. The 5-argument form is present in every VTK
    // version this library supports, so use that one instead of relying
    // on subId being available.
    double p[3] = {point[0], point[1], point[2]};
    vtkGenericCell* genCell = cell ? cell : this->scratchCell;
    return this->locator->FindCell(p, tol2, genCell, pcoords, weights);
}


void
vmtXYZCellLocator::setPeriodicityLengthX(double periodX) {
    if (periodX > 0) {
        std::stringstream msg;
        msg << "vmtXYZCellLocator does not support periodicity (periodX=" << periodX
            << " > 0 requested) -- use vmtLonLatCellLocator for a periodic-longitude grid";
        mntlog::error(__FILE__, __func__, __LINE__, msg.str());
    }
}


void
vmtXYZCellLocator::setCubedSphere(bool isCubedSphere) {
    if (isCubedSphere) {
        mntlog::error(__FILE__, __func__, __LINE__,
            "vmtXYZCellLocator does not support the cubed-sphere spherical-patch model "
            "-- use vmtLonLatCellLocator");
    }
}


bool
vmtXYZCellLocator::invertBilinearPatch(const Vec3& target, const Vec3 verts[4],
                                        double& xi, double& eta) const {
    const int maxIter = 30;
    const double newtonTol2 = 1.e-28; // (1.e-14)^2

    double xi0 = xi, eta0 = eta;
    Vec3 b = this->bilinearMap(xi, eta, verts);
    double tt = dot(target, target);
    if (tt <= 0.) {
        this->closestPointOnBilinearPatch(target, verts, xi, eta);
        return true;
    }
    double lambda = dot(b, target) / tt;

    for (int iter = 0; iter < maxIter; ++iter) {

        // residual and Jacobian columns of F(xi, eta, lambda) = B(xi, eta) - lambda*target
        Vec3 res = this->bilinearMap(xi, eta, verts) - lambda * target;
        Vec3 bXsi = (1. - eta) * (verts[1] - verts[0]) + eta * (verts[2] - verts[3]);
        Vec3 bEta = (1. - xi) * (verts[3] - verts[0]) + xi * (verts[2] - verts[1]);
        Vec3 c = -1.0 * target;

        // solve [bXsi bEta c] (dXsi, dEta, dLambda) = -res by Cramer's rule
        Vec3 bEtaXc = crossProduct(bEta, c);
        double det = dot(bXsi, bEtaXc);
        if (std::abs(det) < 1.e-14 * std::sqrt(dot(bXsi, bXsi) * dot(bEta, bEta) * tt)) {
            // ray (nearly) parallel to the patch
            xi = xi0;
            eta = eta0;
            this->closestPointOnBilinearPatch(target, verts, xi, eta);
            return true;
        }
        Vec3 rhs = -1.0 * res;
        double dXsi = dot(rhs, bEtaXc) / det;
        double dEta = dot(bXsi, crossProduct(rhs, c)) / det;
        double dLambda = dot(bXsi, crossProduct(bEta, rhs)) / det;

        xi += dXsi;
        eta += dEta;
        lambda += dLambda;

        if (dXsi * dXsi + dEta * dEta < newtonTol2) {
            return lambda > 0.;
        }
    }
    return false;
}


bool
vmtXYZCellLocator::closestPointOnBilinearPatch(const Vec3& target, const Vec3 verts[4],
                                                double& xi, double& eta) const {
    const int maxIter = 30;
    const double newtonTol2 = 1.e-24; // (1.e-12)^2
    const double h = 1.e-6;

    for (int iter = 0; iter < maxIter; ++iter) {

        Vec3 res = this->bilinearMap(xi, eta, verts) - target;

        Vec3 rXsi = (this->bilinearMap(xi + h, eta, verts) -
                     this->bilinearMap(xi - h, eta, verts)) / (2. * h);
        Vec3 rEta = (this->bilinearMap(xi, eta + h, verts) -
                     this->bilinearMap(xi, eta - h, verts)) / (2. * h);

        double a11 = dot(rXsi, rXsi), a12 = dot(rXsi, rEta), a22 = dot(rEta, rEta);
        double b1 = -dot(rXsi, res), b2 = -dot(rEta, res);
        double det = a11 * a22 - a12 * a12;
        if (std::abs(det) < 1.e-30) {
            return false;
        }
        double dXsi = (a22 * b1 - a12 * b2) / det;
        double dEta = (a11 * b2 - a12 * b1) / det;

        xi += dXsi;
        eta += dEta;

        if (dXsi * dXsi + dEta * dEta < newtonTol2) {
            return true;
        }
    }
    return false;
}


void
vmtXYZCellLocator::getCellVerts(vtkIdType cellId, Vec3 verts[4]) const {
    vtkCell* c = this->grid->GetCell(cellId);
    vtkPoints* pts = c->GetPoints();
    for (int i = 0; i < 4; ++i) {
        verts[i] = Vec3(pts->GetPoint(i)) - this->projCentre;
    }
}


void
vmtXYZCellLocator::getParametricCoords(vtkUnstructuredGrid* grid, vtkIdType cellId,
                                       const Vec3& p, Vec3& pcoords) {
    Vec3 verts[4];
    this->getCellVerts(cellId, verts);
    double xi = 0.5, eta = 0.5;
    this->invertBilinearPatch(p - this->projCentre, verts, xi, eta);
    pcoords[0] = xi;
    pcoords[1] = eta;
    pcoords[2] = 0.0;
}


void
vmtXYZCellLocator::getEdgePlanes(vtkIdType cellId, Vec3 normals[4]) const {
    Vec3 verts[4];
    this->getCellVerts(cellId, verts);
    Vec3 centroid = 0.25 * (verts[0] + verts[1] + verts[2] + verts[3]);
    for (int k = 0; k < 4; ++k) {
        Vec3 n = crossProduct(verts[k], verts[(k + 1) % 4]);
        double nn = std::sqrt(dot(n, n));
        if (nn > 0.) {
            n = (dot(n, centroid) >= 0. ? 1.0 : -1.0) / nn * n;
        }
        // a degenerate edge (coincident corners) gives n = 0: no constraint
        normals[k] = n;
    }
}


bool
vmtXYZCellLocator::coneContains(const Vec3 normals[4], const Vec3& p, double eps) const {
    double tol = eps * std::sqrt(dot(p, p));
    for (int k = 0; k < 4; ++k) {
        if (dot(normals[k], p) < -tol) {
            return false;
        }
    }
    return true;
}


double
vmtXYZCellLocator::coneExit(const Vec3 normals[4], const Vec3& pBeg, const Vec3& direction,
                            double t0) const {
    // p(t) = pBeg + t*direction is inside plane k while dot(n_k, p(t)) >= 0,
    // which is linear in t: leaves it at t = -dot(n_k, pBeg)/dot(n_k, direction)
    // if dot(n_k, direction) < 0. A line lying in the plane never leaves it.
    const double tiny = 1.e-14 * std::sqrt(dot(direction, direction));
    double tExit = 1.0;
    for (int k = 0; k < 4; ++k) {
        double nd = dot(normals[k], direction);
        if (nd < -tiny) {
            double tk = -dot(normals[k], pBeg) / nd;
            tExit = std::min(tExit, std::max(tk, t0));
        }
    }
    return tExit;
}


vtkIdType
vmtXYZCellLocator::findCellAlongLine(const Vec3& pBeg, const Vec3& direction, double t0,
                                     double& tExit) {
    const double eps = 1.e-12;
    Vec3 p = pBeg + t0 * direction;
    // the cells' bounding boxes hug the surface, while a chord sags below
    // it: search around the radial projection of p onto the sphere through
    // the chord's end points instead
    double r = 0.5 * (std::sqrt(dot(pBeg, pBeg)) + std::sqrt(dot(pBeg + direction, pBeg + direction)));
    double pn = std::sqrt(dot(p, p));
    Vec3 q = (pn > 0. ? r / pn : 1.0) * p;

    vtkIdList* candidates = vtkIdList::New();
    vtkIdType bestCellId = -1;
    double bestExit = t0;
    double halfWidth = 1.e-3 * std::max(r, 1.e-300);
    q += this->projCentre;   // absolute coordinates, for the bucket search
    for (int grow = 0; grow < 12 && bestCellId < 0; ++grow) {
        double bbox[6] = {q[0] - halfWidth, q[0] + halfWidth,
                          q[1] - halfWidth, q[1] + halfWidth,
                          q[2] - halfWidth, q[2] + halfWidth};
        candidates->Reset();
        this->locator->FindCellsWithinBounds(bbox, candidates);
        // of the cells containing p, take the one the line stays in the
        // longest: at a vertex or edge several cells contain p, but only
        // the ones the line actually enters have tExit > t0
        Vec3 normals[4];
        for (vtkIdType i = 0; i < candidates->GetNumberOfIds(); ++i) {
            vtkIdType cId = candidates->GetId(i);
            this->getEdgePlanes(cId, normals);
            if (!this->coneContains(normals, p, eps)) {
                continue;
            }
            double te = this->coneExit(normals, pBeg, direction, t0);
            if (te > bestExit) {
                bestExit = te;
                bestCellId = cId;
            }
        }
        halfWidth *= 4.;
    }
    candidates->Delete();
    tExit = bestExit;
    return bestCellId;
}


std::vector< std::pair<vtkIdType, Vec4> >
vmtXYZCellLocator::findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd) {

    std::vector< std::pair<vtkIdType, Vec4> > result;

    Vec3 direction = pEnd - pBeg;
    const double eps = 10. * std::numeric_limits<double>::epsilon();
    if (dot(direction, direction) < eps) {
        // zero-length segment
        return result;
    }
    // relative to the projection centre
    Vec3 pBegC = pBeg - this->projCentre;

    // Each cell is treated as the cone of rays from the projection centre
    // (see setProjectionCentre) through its
    // bilinear patch (central projection, see invertBilinearPatch). Its
    // boundary is the 4 planes through the centre and the cell's straight
    // edges, so a point is in the cell iff it is on the inner side of all 4
    // planes, and the line pBeg->pEnd leaves the cell where it crosses one
    // of them -- a linear equation, solved exactly. The neighbour across
    // that edge shares the plane, so the line's exit point from one cell is
    // exactly its entry point into the next, and line integrals telescope.
    // This also works for a chord between two points of a sphere, which
    // dips below the cells' patches (a ray-vs-patch test would miss it).

    // safety cap against an infinite loop on degenerate geometry
    const int maxSegments = 100000;
    int numSegments = 0;
    double t0 = 0.0;
    while (t0 < 1.0 - eps) {
        if (++numSegments > maxSegments) {
            std::stringstream msg;
            msg << "exceeded " << maxSegments << " segments tracing a line -- "
                   "likely degenerate geometry";
            mntlog::error(__FILE__, __func__, __LINE__, msg.str());
            break;
        }
        double tExit;
        vtkIdType cellId = this->findCellAlongLine(pBegC, direction, t0, tExit);
        if (cellId < 0) {
            // the line has left the grid
            break;
        }
        Vec4 lambdaInOutPeriodFold;
        lambdaInOutPeriodFold[0] = t0;
        lambdaInOutPeriodFold[1] = tExit;
        lambdaInOutPeriodFold[2] = 0.0; // no periodicity for a plain embedded (x,y,z) mesh
        lambdaInOutPeriodFold[3] = 0.0; // no pole-folding either
        result.push_back(std::pair<vtkIdType, Vec4>(cellId, lambdaInOutPeriodFold));
        t0 = tExit;
    }

    return result;
}
