#include <mntLineGridIntersector.h>
#include <vtkGenericCell.h>

LineGridIntersector::LineGridIntersector(vtkUnstructuredGrid* grid) {

    this->tol = 10 * std::numeric_limits<double>::epsilon();

    this->locator = vtkOBBTree::New();
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
    Vector<double> u = this->direction;
    u *= (1./std::sqrt(lengthSqr));

    // extend the line by a little bit to be sure that the start/end points 
    // are included
    pBeg = pA - 100 * this->tol * u;
    pEnd = pB + 100 * this->tol * u;

    this->tValues.resize(0);

    // add the start point if it is in a cell
    if (this->locator->FindCell( (double*) pa ) >= 0) { // SHOULD WE USE THE VERSION WITH TOL2?
        this->tValues.push_back(0.0);
    }

    //
    // collect the intersection points
    //

    int found = 1;
    double tLocal, tVal;
    double pcoords[3];
    int subId;
    vtkIdType cellId;
    vtkGenericCell* cell = vtkGenericCell::New();
    while (found) {
        found = this->locator->IntersectWithLine(&pBeg[0], &pEnd[0], this->tol, 
                                                 tLocal, &xPoint[0], pcoords, subId, cellId, cell);  // output
        if (found > 0) {
            // compute the linear parameter using the intersection point
            tVal = dot(xPoint - pA, direction)/lengthSqr;
            this->tValues.push_back(tVal);
            // slide the starting point 
            pBeg = xPoint + this->tol * direction;
        }
    }

    // add the end point if it is in a cell
    if (this->locator->FindCell( (double*) pb ) >= 0) { // SHOULD WE USE THE VERSION WITH TOL2?
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
