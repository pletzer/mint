#include "mntPolysegmentIter3d.h"
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <MvVector.h>
#include <limits>

#define DEBUG_PRINT 0


PolysegmentIter3d::PolysegmentIter3d(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                                     const double pa[], const double pb[]) {

    // small tolerances 
    this->eps = 10 * std::numeric_limits<double>::epsilon();
    this->eps100 = 100 * this->eps;

    // store the grid and the grid locator
    this->grid = grid;
    this->locator = locator;

    // reset
    this->cellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);

    // store the begin/end positions
    Vector<double> pA(3), pB(3), direction(3);
    double lengthSqr = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        pA[i] = pa[i];
        pB[i] = pb[i];
        double dist = pb[i] - pa[i];
        direction[i] = dist;
        lengthSqr += dist * dist;
    }

    // length must be > 0
    if (lengthSqr == 0) return;

    double weights[8]; // hex cell
    std::vector<double> xts; // linear parametric coord of the intersection points
    double pcoords[3], pcoords0[3], pcoords1[3]; // parametric coordinates of a point in a cell
    Vector<double> xpoint(3); // intersection point

    // always add the starting point
    xts.push_back(0.);

    // collect the intersection points
    int found = 1;
    double t, tGlobal;
    int subId;
    vtkIdType cellId;
    vtkGenericCell* cell = vtkGenericCell::New();

    Vector<double> pBeg = pA;
    while (found && dot(direction, pB - pBeg) > 0.0) {
        std::cerr << "*** 4. pBeg = " << pBeg << " pB = " << pB << " eps = " << this->eps << " \n";
        found = this->locator->IntersectWithLine(&pBeg[0], &pB[0], this->eps, 
                                                 t, &xpoint[0], pcoords, subId, cellId, cell);  // output
        std::cerr << "*** 4.2 found = " << found << " t = " << t << " xpoint = " << xpoint << "\n";
        if (found > 0) {
            // compute the linear parameter using the intersection point
            tGlobal = dot(xpoint - pA, direction)/lengthSqr;
            std::cerr << "--- tGlobal = " << tGlobal << '\n';
            xts.push_back(tGlobal);

            // slide the starting point 
            pBeg = xpoint + this->eps100 * direction;
        }
    }

    // always add the end point
    xts.push_back(1.0);
    std::cerr << "*** 5 xts = "; for (size_t i = 0; i < xts.size(); ++i) std::cerr << xts[i] << ','; std::cerr << "\n";

    // find all the cells between the t values
    for (size_t iSeg = 0; iSeg < xts.size() - 1; ++iSeg) {

        if (std::abs(xts[iSeg + 1] - xts[iSeg]) < this->eps) {
            continue;
        }

        // make the sure the target points are slightly inside the cell
        Vector<double> xpoint0 = pA + (xts[iSeg + 0] + this->eps100)*direction;
        Vector<double> xpoint1 = pA + (xts[iSeg + 1] - this->eps100)*direction;

        // find the cell index for xpoint0 and xpoint1, should be the same
        vtkIdType cellId0 = this->locator->FindCell(&xpoint0[0], this->eps, cell, pcoords0, weights);
        vtkIdType cellId1 = this->locator->FindCell(&xpoint1[0], this->eps, cell, pcoords1, weights);

        std::cerr << "... cellId0 = " << cellId0 << " pcoords0 = ";
        for (size_t i = 0; i < 3; ++i) std::cerr << pcoords0[i] << ','; 
        std::cerr << " cellId1 = " << cellId1 << " pcoords1 = ";
        for (size_t i = 0; i < 3; ++i) std::cerr << pcoords1[i] << ','; 
        std::cerr << '\n';

        if ( cellId0 == cellId1 ) {
            if (cellId0 >= 0) {
                this->cellIds.push_back(cellId0);
                this->segTas.push_back(xts[iSeg + 0]);
                this->segTbs.push_back(xts[iSeg + 1]);
                this->segXias.push_back( std::vector<double>(pcoords0, pcoords0 + 3) );
                this->segXibs.push_back( std::vector<double>(pcoords1, pcoords1 + 3) );
            }
            else {
                std::cerr << "Warning: cellId0 is not valid " << cellId0 << '\n';
            }
        }
        else {
            std::cerr << "Warning: cellId0 = " << cellId0 << " != cellId1 = " << cellId1 << '\n';
        }
    }

    this->numSegs = this->cellIds.size(); 

    // reset the iterator
    this->reset();

    cell->Delete();
}


double 
PolysegmentIter3d::getIntegratedParamCoord() const {
    if (this->numSegs > 0) {
        return this->segTbs[this->numSegs - 1];
    }
    else {
        return 0;
    }
}


void 
PolysegmentIter3d::reset() {
    this->index = 0;
}


bool
PolysegmentIter3d::next() {
    if (this->index < this->numSegs - 1) {
        this->index++;
        return true;
    }
    return false;
}


vtkIdType 
PolysegmentIter3d::getCellId() const {
    return this->cellIds[this->index];
}

const std::vector<double>& 
PolysegmentIter3d::getBegCellParamCoord() const {
    return this->segXias[this->index];
}        


const std::vector<double>& 
PolysegmentIter3d::getEndCellParamCoord() const {
    return this->segXibs[this->index];
}
 

double
PolysegmentIter3d::getBegLineParamCoord() const {
    return this->segTas[this->index];
}
        

double 
PolysegmentIter3d::getEndLineParamCoord() const {
    return this->segTbs[this->index];
}
     

double 
PolysegmentIter3d::getCoefficient() const {
    return 1;
}
 
size_t
PolysegmentIter3d::getNumberOfSegments() const {
    return this->numSegs;
}


