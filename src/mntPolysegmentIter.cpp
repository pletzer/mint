#include <mntPolysegmentIter.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <limits>
#include <map>
#include <vector>

struct TCmpFunctor {
    TCmpFunctor(const std::vector<double>& ts) {
        this->tVals = ts;
    }
    bool operator()(size_t i, size_t j) {
        return (this->tVals[i] < this->tVals[j]);
    }
    std::vector<double> tVals;
};


PolysegmentIter::PolysegmentIter(vtkUnstructuredGrid* grid,
                                 vmtCellLocator* locator, 
                                 const double p0In[], const double p1In[],
                                 double periodX) {

    // small tolerances 
    this->eps = 10 * std::numeric_limits<double>::epsilon();
    this->eps100 = 100. * this->eps;

    this->periodX = periodX;

    // set the grid and the grid locator
    this->grid = grid;
    this->locator = locator;

    Vec3 p0(p0In);
    this->__makePeriodic(p0);
    Vec3 p1(p1In);
    this->__makePeriodic(p1);

    Vec3 dp = p1 - p0;

    std::vector< std::pair<vtkIdType, Vec3> > cellIdLambdasPeriod = this->locator->findIntersectionsWithLine(p0, p1);

    // arrays of cell Ids, start/end "t" values, start/end "xi" param coords, and 
    // duplicity coefficients for each subsegment
    this->segCellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);
    this->segCoeffs.resize(0);

    double closestPoint[3];
    int subId;
    Vec3 xia, xib;
    double dist2;
    double weights[8];
    for (const auto& cIdLamP : cellIdLambdasPeriod) {

        vtkIdType cId = cIdLamP.first;

        double ta = cIdLamP.second[0];
        double tb = cIdLamP.second[1];
        double periodOffset = cIdLamP.second[2];

        Vec3 pa = p0 + dp*ta;
        Vec3 pb = p0 + dp*tb;
        pa[0] += periodOffset;
        pb[0] += periodOffset;

        // compute the cell parametric coords
        vtkCell* cell = this->grid->GetCell(cId);
        cell->EvaluatePosition(&pa[0], closestPoint, subId, &xia[0], dist2, weights);
        cell->EvaluatePosition(&pb[0], closestPoint, subId, &xib[0], dist2, weights);

        // fill in
        this->segCellIds.push_back(cId);
        this->segTas.push_back(ta);
        this->segTbs.push_back(tb);
        this->segXias.push_back(xia);
        this->segXibs.push_back(xib);
        this->segCoeffs.push_back(1.0);
    }

    // assign coefficients that account for duplicity, ie segments 
    // that are shared between two cells. Output is this->segCoeffs
    this->__assignCoefficientsToSegments();

    this->numSegs = this->segCellIds.size();

    // compute the total, integrated linear param coord
    // should amoount to 1 if the target is entirely 
    // contained within the source grid
    this->totalT = 0.0;
    for (size_t i = 0; i < this->numSegs; ++i) {
        double ta = this->segTas[i];
        double tb = this->segTbs[i];
        double coeff = this->segCoeffs[i];
        this->totalT += (tb - ta) * coeff;
    }

    // reset the iterator
    this->reset();
}


double 
PolysegmentIter::getIntegratedParamCoord() const {
    return this->totalT;
}


void 
PolysegmentIter::reset() {
    this->index = 0;
}


bool
PolysegmentIter::next() {
    if (this->index < this->numSegs - 1) {
        this->index++;
        return true;
    }
    return false;
}


vtkIdType 
PolysegmentIter::getCellId() const {
    return this->segCellIds[this->index];
}

const Vec3& 
PolysegmentIter::getBegCellParamCoord() const {
    return this->segXias[this->index];
}        


const Vec3& 
PolysegmentIter::getEndCellParamCoord() const {
    return this->segXibs[this->index];
}
 

double
PolysegmentIter::getBegLineParamCoord() const {
    return this->segTas[this->index];
}
        

double 
PolysegmentIter::getEndLineParamCoord() const {
    return this->segTbs[this->index];
}
     

double 
PolysegmentIter::getCoefficient() const {
    return this->segCoeffs[this->index];
}
 
size_t
PolysegmentIter::getNumberOfSegments() const {
    return this->numSegs;
}


///////////////////////////////////////////////////////////////////////////////
// private methods

void
PolysegmentIter::__assignCoefficientsToSegments() {

    size_t n = this->segCellIds.size();

    // copy
    std::vector<vtkIdType> sCellIds = this->segCellIds;
    std::vector<double> sTas = this->segTas;
    std::vector<double> sTbs = this->segTbs;
    std::vector<Vec3> sXias = this->segXias;
    std::vector<Vec3> sXibs = this->segXibs;
    std::vector<double> sCoeffs = this->segCoeffs;

    this->segCellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);
    this->segCoeffs.resize(0);
    for (size_t i = 0; i < n; ++i) {
        double ta = sTas[i];
        double tb = sTbs[i];
        if (std::abs(tb - ta) > this->eps100) {
            this->segCellIds.push_back(sCellIds[i]);
            this->segTas.push_back(sTas[i]);
            this->segTbs.push_back(sTbs[i]);
            this->segXias.push_back(sXias[i]);
            this->segXibs.push_back(sXibs[i]);
            this->segCoeffs.push_back(sCoeffs[i]);
        }
    }

    // reduce contribution for the overlapping segments. If two 
    // segments overlap then the coefficient of the first segment
    // is set to 1.0 - overlap/(tb - ta). Assumes the overlap 
    // can only happen for pairs of segment
    n = this->segCellIds.size(); // changed after removing zero length sub-segments
    // iterate over sub-segment pairs
    for (int i0 = 0; i0 < (int) n - 1; ++i0) {
        size_t i1 = i0 + 1;
        double ta0 = this->segTas[i0];
        double tb0 = this->segTbs[i0];
        double ta1 = this->segTas[i1];
        double tb1 = this->segTbs[i1];
        double overlap = std::max(0., std::min(tb0, tb1) - std::max(ta1, ta0));
        this->segCoeffs[i0] = 1.0 - overlap/(tb0 - ta0);
    }

}

void  
PolysegmentIter::__makePeriodic(Vec3& v) {

    // fix start/end points if they fall outside the domain and the domain is periodic
    if (this->periodX > 0.) {
        double xmin = this->grid->GetBounds()[0];
        double xmax = this->grid->GetBounds()[1];
        if (v[0] < xmin) {
            std::cerr << "Warning: adding x periodicity length " << this->periodX << 
                         " to point " << v << "\n";
            v[0] += this->periodX;
        }
        else if (v[0] > xmax) {
            std::cerr << "Warning: subtracting x periodicity length " << this->periodX << 
                         " from point " << v << "\n";
            v[0] -= this->periodX;
        }
    }

}

