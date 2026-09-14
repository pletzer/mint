#include <vmtXYZCellLocator.h>
#include <mntLogger.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <sstream>
#include <limits>
#include <algorithm>


vmtXYZCellLocator::vmtXYZCellLocator() {
    this->locator = vtkStaticCellLocator::New();
    this->scratchCell = vtkGenericCell::New();
    this->grid = nullptr;
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


bool
vmtXYZCellLocator::pointIsInCell(vtkIdType cellId, const double p[3], double tol,
                                  double& xi, double& eta) const {
    Vec3 verts[4];
    vtkCell* c = this->grid->GetCell(cellId);
    vtkPoints* pts = c->GetPoints();
    for (int i = 0; i < 4; ++i) {
        verts[i] = Vec3(pts->GetPoint(i));
    }
    Vec3 target(p);
    // NOTE: deliberately NOT requiring invertBilinearPatch's own strict
    // convergence flag (its Newton iteration stops as soon as the step size
    // is below 1e-12, which is appropriate for a single, near-surface
    // point-location query -- vmtLonLatCellLocator's
    // containsPointCubedSphere, the pattern this was based on, only ever
    // uses it that way). Here, bisection tests points that can start well
    // away from where the warm start is actually valid, and the iteration
    // can legitimately need a few more than maxIter steps, or settle to a
    // slightly looser tolerance, while still landing on the geometrically
    // correct answer -- requiring strict convergence as a PRECONDITION for
    // "inside" was tried first and found, empirically, to produce spurious
    // "outside" verdicts for points that were actually still inside,
    // stalling the walk in scripts/xyz_regrid_prototype.py-sized steps
    // (its own locateInCell never checks convergence either, only the
    // final (xi, eta) bounds -- this now matches that exactly).
    this->invertBilinearPatch(target, verts, xi, eta);
    return xi >= -tol && xi <= 1. + tol && eta >= -tol && eta <= 1. + tol;
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

    // NOTE: this is deliberately NOT a ray-cast (VTK's IntersectWithLine).
    // A straight 3D chord between two points on a curved surface dips
    // strictly BELOW that surface everywhere except at its own endpoints
    // (verified numerically: ~0.0038 units below a unit sphere at the
    // midpoint of a 10-degree-long chord) -- so a literal ray-vs-surface
    // intersection test would find zero intermediate cells for any chord
    // spanning more than one cell width, exactly the "coarser destination
    // edge crossing several finer source cells" case regridding actually
    // needs. Instead: walk along the chord, find the cell containing the
    // current point (via the already-correct, bucket-indexed FindCell),
    // then bisect forward to find where the chord exits that cell's own
    // parametric [0,1]^2 domain, and repeat from there. Prototyped and
    // validated in scripts/xyz_regrid_prototype.py: exact (to numerical
    // noise, ~1e-7) on a flat surface, and converging at roughly
    // O(h^2.5-3.5) on a curved one as resolution increases -- degrading
    // only at deliberately unrealistic resolutions (e.g. a destination
    // cell spanning most of a cubed-sphere panel), which is expected and
    // accepted, not a bug this method needs to work around.
    const double tol = 1.e-6;   // parametric-coordinate ("still inside") tolerance
    // FindCell's own tolerance -- deliberately generous (NOT the tight 1e-6
    // an ordinary point-location query would use). A coarse destination
    // edge's chord can sag well away from the true source surface (e.g.
    // ~0.076 units, for a unit sphere, at the midpoint of a ~45-degree
    // chord -- consistent with the ~0.0038 sag measured for a 10-degree
    // one, since sag grows roughly with the square of the chord's angular
    // length): FindCell would otherwise reject the closest cell outright
    // as "too far" purely because of that sag, not because the point is
    // actually off the grid. The real containment decision here is
    // pointIsInCell's parametric (xi, eta) test below, not physical
    // distance, so this only needs to be generous enough that FindCell's
    // own distance-based rejection never fires first.
    const double tol2 = 1.0;
    const int maxBisections = 60;
    // nudge past a cell boundary once found -- must be well outside the
    // fuzzy "still inside" tolerance zone (tol above) or FindCell can keep
    // re-finding the just-exited cell (bisection alone narrows to near
    // machine precision, so a machine-epsilon-scale nudge does NOT escape
    // that zone -- this was tried first and found empirically to produce
    // nothing but filtered-out zero-length segments). 1000x tol is small
    // relative to any reasonably-sized cell, comfortably larger than tol.
    const double nudge = 1000. * tol;

    // safety cap: a well-formed mesh should need at most a handful of
    // segments per chord; this guards against an infinite loop on
    // pathological/degenerate geometry rather than assuming it can't happen
    const int maxSegments = 100000;

    double t = 0.0;
    int numSegments = 0;
    // persists across iterations -- the last segment's converged (xi, eta),
    // reused as the Gauss-Newton initial guess for the NEXT segment's own
    // cell (matches scripts/xyz_regrid_prototype.py's findIntersectionsWithLineXYZ,
    // which carries xiGuess/etaGuess forward the same way). This matters
    // most for the FindClosestPoint fallback below, which (unlike FindCell)
    // has no parametric coordinates of its own to offer as a guess -- an
    // arbitrary restart at (0.5, 0.5) was tried first and found, empirically,
    // to make invertBilinearPatch diverge for points near a cell's edge or
    // corner (exactly where this fallback tends to fire), truncating the
    // walk well short of t=1.
    double xiPrev = 0.5, etaPrev = 0.5;
    while (t < 1.0 - eps) {

        if (++numSegments > maxSegments) {
            std::stringstream msg;
            msg << "exceeded " << maxSegments << " segments tracing a line -- "
                   "likely degenerate geometry or a bisection tolerance issue";
            mntlog::error(__FILE__, __func__, __LINE__, msg.str());
            break;
        }

        Vec3 p = pBeg + direction * t;
        double pcoords[3], weights[8];
        vtkIdType cellId = this->FindCell(&p[0], tol2, this->scratchCell, pcoords, weights);
        double xi0 = pcoords[0], eta0 = pcoords[1];
        bool foundInside = false;
        if (cellId >= 0) {
            double xiCheck = xi0, etaCheck = eta0;
            if (this->pointIsInCell(cellId, &p[0], tol, xiCheck, etaCheck)) {
                xi0 = xiCheck;
                eta0 = etaCheck;
                foundInside = true;
            }
        }
        if (!foundInside) {
            // Either FindCell's own containment/tol2 test failed outright,
            // or it returned some OTHER, nearby-but-wrong cell (its
            // bucket search can hand back a plausible-looking candidate
            // that our own, more precise (xi, eta) test then rejects --
            // e.g. right after a cell-boundary crossing, where the exact
            // nearest cell isn't necessarily the one whose bucket the
            // query point happens to fall in). Either way this is
            // expected (not "chord left the grid") for a coarse
            // destination edge: the straight chord can sag well below the
            // *source* cells' own tightly-fitting bounding boxes near the
            // surface (each flat quad's bounding box is only as thick as
            // that single cell's own sag, much less than the destination
            // chord's), so the right bucket can easily end up empty, or
            // hand back a neighbour instead.
            //
            // A single "nearest cell" (FindClosestPoint) is not always the
            // RIGHT cell either: right at a boundary crossing, the true owning
            // cell and the just-exited one can be nearly equidistant from
            // p in plain 3D distance (e.g. both bow slightly away from p
            // near a shared corner), so FindClosestPoint can keep handing
            // back the cell we just left. Gather every candidate in a
            // small box around p instead (growing the box until it is
            // non-empty) and test each one with our own (xi, eta) check --
            // this is the bucket-indexed analogue of
            // scripts/xyz_regrid_prototype.py's findContainingCell, which
            // scans every cell rather than trusting a single "closest" one.
            vtkIdList* candidates = vtkIdList::New();
            double boxHalfWidth = 1.e-2;
            const double maxBoxHalfWidth = 2.0; // generous vs. any unit-sphere-scale mesh
            for (int grow = 0; grow < 12 && candidates->GetNumberOfIds() == 0; ++grow) {
                double bbox[6] = {p[0] - boxHalfWidth, p[0] + boxHalfWidth,
                                   p[1] - boxHalfWidth, p[1] + boxHalfWidth,
                                   p[2] - boxHalfWidth, p[2] + boxHalfWidth};
                this->locator->FindCellsWithinBounds(bbox, candidates);
                boxHalfWidth = std::min(boxHalfWidth * 4., maxBoxHalfWidth);
            }
            for (vtkIdType i = 0; i < candidates->GetNumberOfIds() && !foundInside; ++i) {
                vtkIdType candCellId = candidates->GetId(i);
                // try the persisted warm start first, and a neutral
                // center-of-cell restart second -- a point landing very
                // close to a shared CORNER (as opposed to just an edge) of
                // the source mesh can leave invertBilinearPatch's
                // Gauss-Newton in an ill-conditioned spot when warm-started
                // from a value that itself sits right at a cell boundary
                // (e.g. xi or eta == 1 exactly); a plain (0.5, 0.5) restart
                // has no such degeneracy and empirically recovers those
                // cases.
                double xiCheck = xiPrev, etaCheck = etaPrev;
                bool inside2 = this->pointIsInCell(candCellId, &p[0], tol, xiCheck, etaCheck);
                if (!inside2) {
                    xiCheck = 0.5;
                    etaCheck = 0.5;
                    inside2 = this->pointIsInCell(candCellId, &p[0], tol, xiCheck, etaCheck);
                }
                if (inside2) {
                    cellId = candCellId;
                    xi0 = xiCheck;
                    eta0 = etaCheck;
                    foundInside = true;
                }
            }
            candidates->Delete();
        }
        if (!foundInside) {
            // neither FindCell nor the nearest-cell fallback found a cell
            // that actually claims this point in its own parametric
            // bounds -- the chord has genuinely left the grid
            break;
        }
        xiPrev = xi0;
        etaPrev = eta0;

        // FIXED warm start for this whole segment (NOT updated per bisection
        // test) -- the converged answer at t. Matches
        // scripts/xyz_regrid_prototype.py's stillInside exactly: that
        // closure captures (xi, eta) once, from the initial containing-cell
        // search, and reuses the SAME values for every bisection test.
        // Updating the warm start to the latest "confirmed inside" solution
        // during bisection was tried first (a seemingly reasonable idea --
        // warm-start from the closest known point) and found, empirically,
        // to make convergence WORSE near a cell boundary: bisection's mid
        // values don't move monotonically towards the true boundary, so an
        // evolving warm start can end up seeded from a point that is
        // (in xi, eta space) no closer to the current test point than a
        // fixed one, while accumulating whatever numerical drift the
        // updates introduce -- observed as the walk crawling forward in
        // steps barely larger than `nudge`, never reaching the cell's true
        // exit.

        Vec3 pEndOfChord = pBeg + direction * 1.0;
        double tb;
        double xiTry = xi0, etaTry = eta0;
        if (this->pointIsInCell(cellId, &pEndOfChord[0], tol, xiTry, etaTry)) {
            // the rest of the chord, all the way to t=1, stays in this cell
            tb = 1.0;
        }
        else {
            // bisect between t (known inside) and 1 (known outside) to
            // find where the chord exits this cell
            double lo = t, hi = 1.0;
            for (int iter = 0; iter < maxBisections; ++iter) {
                double mid = 0.5 * (lo + hi);
                Vec3 pMid = pBeg + direction * mid;
                xiTry = xi0;
                etaTry = eta0;
                if (this->pointIsInCell(cellId, &pMid[0], tol, xiTry, etaTry)) {
                    lo = mid;
                }
                else {
                    hi = mid;
                }
            }
            tb = lo;
        }
        Vec4 lambdaInOutPeriodFold;
        lambdaInOutPeriodFold[0] = t;
        lambdaInOutPeriodFold[1] = tb;
        lambdaInOutPeriodFold[2] = 0.0; // no periodicity for a plain embedded (x,y,z) mesh
        lambdaInOutPeriodFold[3] = 0.0; // no pole-folding either
        result.push_back(std::pair<vtkIdType, Vec4>(cellId, lambdaInOutPeriodFold));

        // nudge past the boundary before locating the next cell
        t = tb + nudge;
    }

    return result;
}
