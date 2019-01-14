#include <mntLineGridIntersector.h>
#include <vtkGenericCell.h>

LineGridIntersector::LineGridIntersector(vtkUnstructuredGrid* grid) {

    this->tol = 10 * std::numeric_limits<double>::epsilon();

    this->locator = vtkCellLocator::New();
    this->locator->SetDataSet(grid);
    this->locator->BuildLocator();
}

LineGridIntersector::~LineGridIntersector() {
    this->locator->Delete();
}

void 
LineGridIntersector::setLine(const double pa[], const double pb[]) {

    this->pA.resize(3);
    this->direction.resize(3);
    Vector<double> pBeg(3), pEnd(3), xPoint(3), pB(3);
    double lengthSqr = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        this->direction[i] = pb[i] - pa[i];
        lengthSqr += this->direction[i] * this->direction[i];
        this->pA[i] = pa[i];
        pB[i] = pb[i];
    }

    // extend the line by a little bit to be sure that the start/end points 
    // are included
    pBeg = pA - 100 * this->tol * this->direction;
    pEnd = pB + 100 * this->tol * this->direction;

    this->tValues.clear();
    vtkIdType cellId;

    // add the start point if it is in a cell
    cellId = this->locator->FindCell( (double*) pa ); // SHOULD WE USE THE VERSION WITH TOL2?
    if (cellId >= 0) { 
        this->tValues.push_back(0.0);
    }

    if (lengthSqr == 0) {
        // not really a line, just a point
        return;
    }

    //
    // collect the intersection points
    //

    int found = 1;
    double tLocal, tVal;
    double pcoords[3];
    int subId;
    vtkGenericCell* cell = vtkGenericCell::New();
    while (found) {
        found = this->locator->IntersectWithLine(&pBeg[0], &pEnd[0], this->tol, 
                                                 tLocal, &xPoint[0], pcoords, subId, cellId, cell);  // output
        if (found > 0) {
            // compute the linear parameter using the intersection point
            tVal = dot(xPoint - pA, direction)/lengthSqr;
            // add if different from previous
            size_t nValues = this->tValues.size();
            double absDiff = std::abs( tVal - this->tValues[nValues - 1] );
            if (nValues == 0 || absDiff > this->tol) {
                this->tValues.push_back(tVal);
            }
            // slide the starting point 
            pBeg = xPoint + this->tol * this->direction;
        }
    }
    cell->Delete();

    // add the end point if it is in a cell
    cellId = this->locator->FindCell( (double*) pb ); // SHOULD WE USE THE VERSION WITH TOL2?
    double absDiff = std::numeric_limits<double>::max();
    size_t nValues = this->tValues.size();
    if (nValues > 0) {
        absDiff = std::abs( 1.0 - this->tValues[nValues - 1] );
    }
    if (cellId >= 0 && absDiff > this->tol) { 
        this->tValues.push_back(1.0);
    }

}

const std::vector<double>& 
LineGridIntersector::getIntersectionLineParamCoords() const {
    return this->tValues;
}

std::vector< Vector<double> > 
LineGridIntersector::getIntersectionPoints() const {
    std::vector< Vector<double> > res;
    res.reserve(this->tValues.size());
    for (size_t i = 0; i < this->tValues.size(); ++i) {
        Vector<double> xPoint = this->direction;
        xPoint *= this->tValues[i];
        xPoint += this->pA;
        res.push_back(xPoint);
    }
    return res;
}

const Vector<double>& 
LineGridIntersector::getDirection() const {
    return this->direction;
}

const Vector<double>& 
LineGridIntersector::getStartPoint() const {
    return this->pA;
}


