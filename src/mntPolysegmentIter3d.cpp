#include "mntPolysegmentIter3d.h"
#include "mntLineGridIntersector.h"
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <MvVector.h>
#include <limits>

#define DEBUG_PRINT 1


PolysegmentIter3d::PolysegmentIter3d(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                                     const double pa[], const double pb[]) {

    // small tolerances 
    this->eps = 10 * std::numeric_limits<double>::epsilon();

    // store the grid and the grid locator
    this->grid = grid;
    this->locator = locator;

    vtkGenericCell* cell;
    double pcoords0[3];
    double pcoords1[3];
    double weights[8];
    int subId;
    double dist2;
    int inside0, inside1; // 1=inside, 0=outside, -1=error

    // reset
    this->cellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);

    LineGridIntersector intersector(grid);

    intersector.setLine(pa, pb);
    const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();

    if (tValues.size() == 0) {
        // the line does not intersect with the grid
#ifdef DEBUG_PRINT
        std::cerr << "Warning: no intersection\n";
#endif
        return;
    }
    if (tValues.size() == 1) {
        // zero length overlap between the line and the grid
#ifdef DEBUG_PRINT
        std::cerr << "Warning: zero length overlap between line and grid\n";
#endif
        return;
    }
    std::vector< Vector<double> > xpoints = intersector.getIntersectionPoints();

    // find all the cells between the t values
    for (size_t iSeg = 0; iSeg < tValues.size() - 1; ++iSeg) {

        if (std::abs(tValues[iSeg + 1] - tValues[iSeg]) < this->eps) {
            // skip if the segment has zero length in t-space
#ifdef DEBUG_PRINT
            std::cerr << "Warning: duplicate intersection points t = " << tValues[iSeg] << ", " << tValues[iSeg + 1] << "\n";
#endif
            continue;
        }

        const Vector<double>& pA = intersector.getStartPoint();
        const Vector<double>& direction = intersector.getDirection();

        double tDiff = tValues[iSeg + 1] - tValues[iSeg + 0];
        double tMid = 0.5*(tValues[iSeg + 0] + tValues[iSeg + 1]);
        Vector<double> pMid = direction;
        pMid *= tMid;
        pMid += pA;

        vtkIdType cellId = this->locator->FindCell(&pMid[0]);
        if (cellId < 0) {
            std::cerr << "Warning: could not find cell at seg mid point t = " 
                      << tMid << " point = " << pMid << " (cellId = " << cellId << ")\n";
            continue;

        }

        // paramatric coords at the start of the segment
        Vector<double> pt = direction;
        pt *= tValues[iSeg + 0];
        pt += pA;
        inside0 = this->grid->GetCell(cellId)->EvaluatePosition(&pt[0], NULL, subId, pcoords0, dist2, weights);
        if (inside0 != 1) {
            std::cerr << "Warning: could not find pcoords at seg start t = " 
                      << tValues[iSeg + 0] << " point = " << pt << " code = " << inside0 << '\n';
            continue;
        }

        // parametric coords at the end of the segment
        pt += tDiff * direction;
        inside1 = this->grid->GetCell(cellId)->EvaluatePosition(&pt[0], NULL, subId, pcoords1, dist2, weights);
        if (inside1 != 1) {
            std::cerr << "Warning: could not find pcoords at seg end t = " 
                      << tValues[iSeg + 1] << " point = " << pt << " code = " << inside0 << '\n';
            continue;
        }

        this->cellIds.push_back(cellId);
        this->segTas.push_back(tValues[iSeg + 0]);
        this->segTbs.push_back(tValues[iSeg + 1]);
        this->segXias.push_back( std::vector<double>(pcoords0, pcoords0 + 3) );
        this->segXibs.push_back( std::vector<double>(pcoords1, pcoords1 + 3) );

    }

    this->numSegs = this->cellIds.size(); 

    // reset the iterator
    this->reset();

    //cell->Delete();

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


