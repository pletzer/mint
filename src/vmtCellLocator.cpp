#include <vmtCellLocator.h>
#include <vtkQuad.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <mntLineLineIntersector.h>
#include <algorithm>

struct LambdaBegFunctor {
    // compare two elements of the array
    bool operator()(const std::pair<vtkIdType, Vec3>& x, 
                    const std::pair<vtkIdType, Vec3>& y) {
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


bool 
vmtCellLocator::containsPoint(vtkIdType faceId, const double point[3], double tol) const {

    tol = std::abs(tol);
    bool res = true;
    double adeltas0[3];
    double adeltas1[3];
    double* val;
    int k;

    std::vector<Vec3> nodes = this->getFacePoints(faceId);
    size_t npts = nodes.size();

    for (size_t i0 = 0; i0 < npts; ++i0) {

        size_t i1 = (i0 + 1) % npts;

        double* p0 = &nodes[i0][0];
        double* p1 = &nodes[i1][0];

        // vector from point to the vertices
        double dx0 = point[0] - p0[0];
        double dx1 = point[0] - p1[0];
        double dy0 = point[1] - p0[1];
        double dy1 = point[1] - p1[1];

        // add/substract periodicity length to minimize the distance between 
        // longitudes
        size_t numPer = this->modPeriodX.size(); // either 1 or 3
        for (size_t j = 0; j < numPer; ++j) {
            // modPeriodX could be {0, -360, 360}
            adeltas0[j] = std::abs(dx0 + this->modPeriodX[j]);
            adeltas1[j] = std::abs(dx1 + this->modPeriodX[j]);
        }

        val = std::min_element(adeltas0, adeltas0 + numPer);
        k = (int) std::distance(adeltas0, val);
        dx0 += this->modPeriodX[k];

        val = std::min_element(adeltas1, adeltas1 + numPer);
        k = (int) std::distance(adeltas1, val);
        dx1 += this->modPeriodX[k];

        double cross = dx0*dy1 - dy0*dx1;
        if (cross < -tol) {
            // negative area
            res = false;
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
//                    std::cerr << "*** adding cell id " << faceId << 
//                                 " m, mLo, mHi = " << m << ',' << mLo << ',' << mHi << 
//                                 " n, nLo, nHi = " << n << ',' << nLo << ',' << nHi << '\n';
                }
            }
        }
    }

}


std::vector< std::pair<vtkIdType, Vec3> >
vmtCellLocator::findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd) {

    Vec3 direction = pEnd - pBeg;
    double p0[] = {pBeg[0], pBeg[1], pBeg[2]};
    double p1[] = {pEnd[0], pEnd[1], pEnd[2]};

    double lambdaInOutPeriod[3];

    // store result
    std::vector< std::pair<vtkIdType, Vec3> > res;

    const double eps = 10 * std::numeric_limits<double>::epsilon();

    vtkIdList* cellIds = vtkIdList::New();

    for (double modPx : this->modPeriodX) {

        p0[0] = pBeg[0] + modPx;
        p1[0] = pEnd[0] + modPx;

        // collect the cell Ids intersected by the line  
        this->FindCellsAlongLine(p0, p1, eps, cellIds);

        // iterate over the intersected cells
        for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {

           vtkIdType cellId = cellIds->GetId(i);

            std::vector<double> lambdas = this->collectIntersectionPoints(cellId, p0, direction);

            if (lambdas.size() >= 2) {

                lambdaInOutPeriod[0] = lambdas[0];
                lambdaInOutPeriod[1] = lambdas[lambdas.size() - 1];
                lambdaInOutPeriod[2] = modPx;

                // found entry/exit points so add
                res.push_back(  std::pair<vtkIdType, Vec3>( cellId, Vec3(lambdaInOutPeriod) )  );
            }
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


