#include <vmtCellLocator.h>
#include <vtkQuad.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <mntLineLineIntersector.h>
#include <algorithm>


inline bool isPointInQuad(const Vec3& targetPoint, std::vector<Vec3>& nodes, double tol) {

    bool res = true;

    // nuber of points in the quad
    size_t npts = nodes.size();

    // iterate over the edges of the quad
    for (size_t i0 = 0; i0 < npts; ++i0) {

        size_t i1 = (i0 + 1) % npts;

        // starting/end points of the edge
        double* p0 = &nodes[i0][0];
        double* p1 = &nodes[i1][0];

        // vectors from point to the vertices
        double dx0 = p0[0] - targetPoint[0];
        double dx1 = p1[0] - targetPoint[0];
        double dy0 = p0[1] - targetPoint[1];
        double dy1 = p1[1] - targetPoint[1];

        // area is positive if point is inside the quad
        double cross = dx0*dy1 - dy0*dx1;
        res &= (cross > -tol);
    }

    return res;
}


struct LambdaBegFunctor {
    // compare two elements of the array
    bool operator()(const std::pair<vtkIdType, Vec4>& x, 
                    const std::pair<vtkIdType, Vec4>& y) {
        return (x.second[0] < y.second[0]);
    }
};


vmtCellLocator::vmtCellLocator() {

    this->grid = NULL;

    // will be determined in SetDataSet
    this->numBucketsX = 1;
    this->numBucketsY = 1;

    double big = std::numeric_limits<double>::max();
    for (size_t i = 0; i < 3; ++i) {
        this->xmin[i] = big;
        this->xmax[i] = -big;
    }
    this->modPeriodX.push_back(0.0);
    this->kFolding.push_back(0);

    // now folding is turned on by default
    this->enableFolding();

}


void 
vmtCellLocator::SetDataSet(vtkUnstructuredGrid* grid) {

    this->grid = grid;

    double* bounds = grid->GetBounds();
    this->xmin[0] = bounds[0];
    this->xmax[0] = bounds[1];
    this->xmin[1] = bounds[2];
    this->xmax[1] = bounds[3];
    this->xmin[2] = bounds[4];
    this->xmax[2] = bounds[5];

    this->lambdaMid = 0.5*(this->xmin[0] + this->xmax[0]);

    // want the buckets to be larger than the cells
    vtkIdType numCells = grid->GetNumberOfCells();
    if (numCells > 5) {
        this->numBucketsX = 5; 
    }
    this->numBucketsY = std::max(1, static_cast<int>(numCells/(this->numBucketsX * 128)));
}

void 
vmtCellLocator::SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) {

    vtkIdType numCells = this->grid->GetNumberOfCells();
    this->numBucketsY = std::max(1, static_cast<int>(numCells/(this->numBucketsX * avgNumFacesPerBucket)));
}


void 
vmtCellLocator::BuildLocator() {

    std::set<vtkIdType> empty;
    // attach an empty set of face Ids to each bucket
    for (int m = 0; m < this->numBucketsX; ++m) {
        for (int n = 0; n < this->numBucketsY; ++n) {
            int bucketId = m * this->numBucketsY + n;
            // create empty bucket
            this->bucket2Faces.insert( std::pair< int, std::set<vtkIdType> >(bucketId, empty) );
        }
    }

    // assign each face to one or more buckets depending on where the face's nodes fall
    // WARNING: this could fail if the buckets are much smaller than some cells!
    vtkIdType numFaces = this->grid->GetNumberOfCells();
    for (vtkIdType faceId = 0; faceId < numFaces; ++faceId) {
        std::vector<Vec3> nodes = getFacePoints(faceId);
        for (const Vec3& p : nodes) {
            // assumes that the buckets are always bigger than the faces, or alternatively
            // that no corners are fully inside a face
            int bucketId = this->getBucketId(&p[0]);
            this->bucket2Faces[bucketId].insert(faceId);
        }
    }

}


void 
vmtCellLocator::setPeriodicityLengthX(double periodX) {

    this->modPeriodX.resize(0);
    this->modPeriodX.push_back(0.0);
    if (periodX > 0) {
        this->modPeriodX.push_back(-periodX);
        this->modPeriodX.push_back(+periodX);
    }

}

void
vmtCellLocator::enableFolding() {
    this->kFolding.resize(2);
    this->kFolding[0] = 0;
    this->kFolding[1] = 1;
}

bool 
vmtCellLocator::containsPoint(vtkIdType faceId, const double point[3], double tol) const {

    tol = std::abs(tol);
    std::vector<Vec3> nodes = this->getFacePoints(faceId);
    Vec3 targetPoint(point);

    bool res = isPointInQuad(targetPoint, nodes, tol);

    return res;
}


bool 
vmtCellLocator::containsPointMultiValued(vtkIdType faceId, const double point[3], double tol) const {

    bool res = false;
    tol = std::abs(tol);
    std::vector<Vec3> nodes = this->getFacePoints(faceId);
    Vec3 targetPoint(point);

    // store the original values
    double lon = targetPoint[0];
    double lat = targetPoint[1];

    // add/substract periodicity length and apply folding if need be
    for (auto kF : this->kFolding) {

        if (kF > 0) {
            // apply folding
            this->foldAtPole(&targetPoint[0]);
        }

        for (const auto& periodX : this->modPeriodX) {

            // add periodicity length
            targetPoint[0] += periodX;
            // point is inside the quad if res is positive for any
            // periodicity lengths
            res |= isPointInQuad(targetPoint, nodes, tol);

            // back to the original values
            targetPoint[0] -= periodX;
        }
    }

    return res;
}


vtkIdType
vmtCellLocator::FindCell(const double point[3], double tol, vtkGenericCell *notUsed, double pcoords[3], double *weights) {

    int bucketId = this->getBucketId(point);
    double closestPoint[3];
    int subId;
    double dist2;

    const std::set<vtkIdType>& faces = this->bucket2Faces.find(bucketId)->second;

    for (const vtkIdType& cId : faces) {
        if (this->containsPoint(cId, point, tol)) {
            vtkCell* quad = this->grid->GetCell(cId);
            quad->EvaluatePosition((double*) point, closestPoint, subId, pcoords, dist2, weights);
            return cId;
        }
    }

    // failed to find cell
    return -1;

}


void
vmtCellLocator::FindCellsAlongLine(const double p0[3], const double p1[3], double tol2, vtkIdList *cellIds) {

    cellIds->Reset();

    int begM, endM, begN, endN, bucketId, begBucketId, endBucketId;

    Vec3 point0(p0);
    Vec3 point1(p1);

    // choose the number of sections heuristically. Too few and we'll end up adding too many
    // cells. No point in having more sections than the number of buckets
    this->getBucketIndices(this->getBucketId(p0), &begM, &begN);
    this->getBucketIndices(this->getBucketId(p1), &endM, &endN);

    int mLo = std::min(begM, endM);
    int mHi = std::max(begM, endM);
    int nLo = std::min(begN, endN);
    int nHi = std::max(begN, endN);

    // dm and dn are positive
    int dm = mHi - mLo + 1;
    int dn = nHi - nLo + 1;

    // break the line into segments, the number of segments should not affect the result. Only used as a performance 
    // improvement. Want more segments when the line is 45 deg. Want more segments when the points are far apart. 
    size_t nSections = std::max(1, std::min(dn, dm));
    Vec3 du = point1 - point0;
    du /= (double) nSections;

    for (size_t iSection = 0; iSection < nSections; ++iSection) {

        // start/end points of the segment
        Vec3 pBeg = point0 + (double) iSection * du;
        Vec3 pEnd = pBeg + du;
    
        // get the start bucket
        begBucketId = this->getBucketId(&pBeg[0]);
        this->getBucketIndices(begBucketId, &begM, &begN);

        // get end bucket
        endBucketId = this->getBucketId(&pEnd[0]);
        this->getBucketIndices(endBucketId, &endM, &endN);

        mLo = std::min(begM, endM);
        mHi = std::max(begM, endM);
        nLo = std::min(begN, endN);
        nHi = std::max(begN, endN);

        // iterate over the buckets
        for (int m = mLo; m <= mHi; ++m) {
            for (int n = nLo; n <= nHi; ++n) {
                bucketId = m * numBucketsY + n;
                for (const vtkIdType& faceId : this->bucket2Faces.find(bucketId)->second) {
                    cellIds->InsertUniqueId(faceId);
                }
            }
        }
    }

}


std::vector< std::pair<vtkIdType, Vec4> >
vmtCellLocator::findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd) {

    Vec3 direction;

    Vec3 p0 = pBeg;
    Vec3 p1 = pEnd;

    Vec4 lambdaInOutPeriodFold;

    // store result
    std::vector< std::pair<vtkIdType, Vec4> > res;

    const double eps = 10 * std::numeric_limits<double>::epsilon();

    vtkIdList* cellIds = vtkIdList::New();

    double lam0 = pBeg[0];
    double the0 = pBeg[1];
    double lam1 = pEnd[0];
    double the1 = pEnd[1];

    for (int kFold : this->kFolding) {

        if (kFold == 1) {

            if (std::abs(p0[1]) <= 90 && std::abs(p1[1]) <= 90) {
                // no need to fold
                continue;
            }

            // transform
            foldAtPole(&p0[0]);
            foldAtPole(&p1[0]);
        }

        direction = p1 - p0;
        double distSq = dot(direction, direction);
        if (distSq < eps) {
            // zero length, skip
            std::cout << "... zero distance: distSq = " << distSq << '\n';
            continue;
        }

        for (double modPx : this->modPeriodX) {

            // add/subtract periodicity to start/end points
            p0[0] += modPx;
            p1[0] += modPx;

            // collect the cell Ids intersected by the line  
            this->FindCellsAlongLine(&p0[0], &p1[0], eps, cellIds);

            // iterate over the intersected cells
            for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {

                vtkIdType cellId = cellIds->GetId(i);

                std::vector<double> lambdas = this->collectIntersectionPoints(cellId, p0, direction);
                size_t nLam = lambdas.size();
                double lamBeg = lambdas.front();
                double lamEnd = lambdas.back();

                if (nLam >= 2 && std::abs(lamEnd - lamBeg) > eps) {

                    lambdaInOutPeriodFold[0] = lamBeg;
                    lambdaInOutPeriodFold[1] = lamEnd;
                    lambdaInOutPeriodFold[2] = modPx;
                    lambdaInOutPeriodFold[3] = kFold;

                    // found entry/exit points so add
                    res.push_back(  std::pair<vtkIdType, Vec4>(cellId, lambdaInOutPeriodFold)  );
                }
            }

            // revert to the original values
            p0[0] = lam0; p0[1] = the0;
            p1[0] = lam1; p1[1] = the1;
        }
    }

    cellIds->Delete();

    // sort by starting lambda
    std::sort(res.begin(), res.end(), LambdaBegFunctor());

    // to avoid double counting, shift lambda entry to be always >= to 
    // the preceding lambda exit and make sure lambda exit >= lambda entry
    for (size_t i = 1; i < res.size(); ++i) {

        double thisLambdaBeg = res[i].second[0];
        double thisLambdaEnd = res[i].second[1];
        double precedingLambdaEnd = res[i - 1].second[1];

        thisLambdaBeg = std::min(thisLambdaEnd, std::max(thisLambdaBeg, precedingLambdaEnd));

        // reset lambda entry
        res[i].second[0] = thisLambdaBeg;
    }

    return res;
}


void 
vmtCellLocator::printBuckets() const {
    for (const auto& b2f : this->bucket2Faces) {
        int bucketId = b2f.first;
        int m, n;
        this->getBucketIndices(bucketId, &m, &n);
        std::cout << "bucket " << bucketId << " (" << m << ',' << n << ") contains faces ";
        for (const auto& faceId : b2f.second) {
            std::cout << faceId << ' ';
        }
        std::cout << '\n';
    }
}


std::vector<double>
vmtCellLocator::collectIntersectionPoints(vtkIdType cellId, 
                                          const Vec3& pBeg,
                                          const Vec3& direction) {

    Vec3 pEnd = pBeg + direction;

    std::vector<double> lambdas;
    // expect two values
    lambdas.reserve(2);

    const double eps = 10 * std::numeric_limits<double>::epsilon();
    const double eps100 = 100*eps;


    std::vector<Vec3> nodes = this->getFacePoints(cellId);
    size_t npts = nodes.size();

    // computes the intersection point of two lines
    LineLineIntersector intersector;

    // is the starting point inside the cell?
    if (this->containsPoint(cellId, &pBeg[0], eps)) {
        lambdas.push_back(0.);
    }

    // iterate over the cell's edges
    size_t i0, i1;
    for (i0 = 0; i0 < npts; ++i0) {

        i1 = (i0 + 1) % npts;

        // compute the intersection point
        double* p0 = &nodes[i0][0];
        double* p1 = &nodes[i1][0];

        intersector.setPoints(2, &pBeg[0], &pEnd[0], &p0[0], &p1[0]);

        if (! intersector.hasSolution(eps)) {
            // no solution, skip
            continue;
        }

        // we have a solution but it could be degenerate

        if (std::abs(intersector.getDet()) > eps) {
            // normal intersection, 1 solution
            Vec2 sol = intersector.getSolution();
            double lambRay = sol[0];
            double lambEdg = sol[1];

            // is it valid? Intersection must be within (p0, p1) and (q0, q1)
            if (lambRay >= (0. - eps100) && lambRay <= (1. + eps100)  && 
                lambEdg >= (0. - eps100) && lambEdg <= (1. + eps100)) {
                // add to list
                lambdas.push_back(lambRay);
            }
        }
        else {
            // det is almost zero
            // looks like the two lines (p0, p1) and (q0, q1) are overlapping
            // add the starting/end points
            const std::pair<double, double>& sol = intersector.getBegEndParamCoords();
            // add start/end linear param coord along line
            lambdas.push_back(sol.first);
            lambdas.push_back(sol.second);

        }
    }

    // is the end point inside the cell?
    if (this->containsPoint(cellId, &pEnd[0], eps)) {
        lambdas.push_back(1.);
    }

    // order the lambdas
    std::sort(lambdas.begin(), lambdas.end());

    return lambdas;
}


