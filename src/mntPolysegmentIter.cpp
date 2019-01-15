#include <mntPolysegmentIter.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <MvVector.h>
#include <limits>

struct TCmpFunctor {
    TCmpFunctor(const std::vector<double>& ts) {
        this->tVals = ts;
    }
    bool operator()(size_t i, size_t j) {
        return (this->tVals[i] < this->tVals[j]);
    }
    std::vector<double> tVals;
};


PolysegmentIter::PolysegmentIter(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                                 const double p0[], const double p1[]) {

    // small tolerances 
    this->eps = 10 * std::numeric_limits<double>::epsilon();
    this->eps100 = 100. * this->eps;
    this->tol = 1.e-3; // to determine if a point is inside a cell

    // set the grid and the grid locator
    this->grid = grid;
    this->locator = locator;

    // cellIds, xis and ts are output
    this->cellIds.resize(0); // cell of each intersection point
    this->xis.resize(0);     // cell parametric coords for each intersection point
    this->ts.resize(0);      // linear param coord for each intersction point
    this->__collectLineGridSegments(p0, p1);

    // gather the intersection points attached to a cell: cellId -> [indx0, indx1, ...] 
    // indx is index in the cellIds, xis and ts arrays
    std::map< vtkIdType, std::vector<size_t> > c2Inds;
    for (size_t i = 0; i < this->cellIds.size(); ++i) {
        vtkIdType cId = this->cellIds[i];
        std:map< vtkIdType, std::vector<size_t> >::iterator it = c2Inds.find(cId);
        if (it != c2Inds.end()) {
            // push_back into existing list value
            it->second.push_back(i);
        }
        else {
            // create new list value entry
            std::vector<size_t> index(1, i);
            std::pair< vtkIdType, std::vector<size_t> > p(cId, index);
            c2Inds.insert(p);
        }
    }

    //
    // build the subsegments
    //

    // arrays of cell Ids, start/end t values, start/end xi param coords, and 
    // duplicity coefficients for each subsegment
    this->segCellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);
    this->segCoeffs.resize(0);

    // iterate over all the cells for which we have intersection points
    for (std::map< vtkIdType, std::vector<size_t> >::const_iterator it = c2Inds.begin();
        it != c2Inds.end(); ++it) {

        // cell Id
        vtkIdType cId = it->first;

        // index array into this->ts, this->cellIds and this->xis
        const std::vector<size_t>& inds = it->second;

        // compute the corresponding linear param coord t's
        size_t n = inds.size();
        std::vector<size_t> iVals(n);
        std::vector<double> tVals(n);
        for (size_t i = 0; i < n; ++i) {
            iVals[i] = i;
            tVals[i] = this->ts[ inds[i] ];
        }

        // sort the indices inds by t values
        std::sort(iVals.begin(), iVals.end(), TCmpFunctor(tVals));

        // create subsegments. Each subsegment has start and end points. Both
        // the start/end points are in the same cell.
        for (int j = 0; j < (int) n - 1; ++j) {
            // indices into the inds lists
            size_t i0 = iVals[j    ];
            size_t i1 = iVals[j + 1];
            // indices into this->ts, this->xis... 
            size_t ia = inds[i0];
            size_t ib = inds[i1];
            double ta = this->ts[ia];
            double tb = this->ts[ib];
            const std::vector<double>& xia = this->xis[ia];
            const std::vector<double>& xib = this->xis[ib];

            // add the cell index, start linear param coord, etc.
            this->segCellIds.push_back(cId);
            this->segTas.push_back(ta);
            this->segTbs.push_back(tb);
            this->segXias.push_back(xia);
            this->segXibs.push_back(xib);
            // will deal with duplicity later
            this->segCoeffs.push_back(1.0);
        }
    }

    // sort all the segments by start linear param coord t values

    size_t n = this->segCellIds.size();
    std::vector<size_t> iVals(n);
    std::vector<double> tVals(n);
    for (size_t i = 0; i < n; ++i) {
        iVals[i] = i;
        tVals[i] = this->segTas[i];
    }

    std::sort(iVals.begin(), iVals.end(), TCmpFunctor(tVals));

    // copy
    std::vector<vtkIdType> sCellIds = this->segCellIds;
    std::vector<double> sTas = this->segTas;
    std::vector<double> sTbs = this->segTbs;
    std::vector< std::vector<double> > sXias = this->segXias;
    std::vector< std::vector<double> > sXibs = this->segXibs;

    // no need to sort this->segCoeffs since all the values are one
    for (size_t j = 0; j < n; ++j) {
        size_t i = iVals[j];
        this->segCellIds[j] = sCellIds[i];
        this->segTas[j] = sTas[i];
        this->segTbs[j] = sTbs[i];
        this->segXias[j] = sXias[i];
        this->segXibs[j] = sXibs[i];
    }

    // assign coefficients that account for duplicity, ie segments 
    // that are shared between two cells. Output is this->segCoeffs
    this->__assignCoefficientsToSegments();

    this->numSegs = this->segCellIds.size();

    // compute the total, integrated linear param coord
    // should amoount to 1 is the target is entirely 
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

const std::vector<double>& 
PolysegmentIter::getBegCellParamCoord() const {
    return this->segXias[this->index];
}        


const std::vector<double>& 
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
    std::vector< std::vector<double> > sXias = this->segXias;
    std::vector< std::vector<double> > sXibs = this->segXibs;
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

    // reduce contribution for overlapping segments. If two 
    // segments overlap then the coefficient of first segment
    // is set to 1.0 - overlap/(tb - ta). Assumes overlap 
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
PolysegmentIter::__collectIntersectionPoints(const double pBeg[], 
                                             const double pEnd[],
                                             std::vector<vtkIdType>& cIds,
                                             std::vector<double>& lambRays,
                                             std::vector< std::vector<double> >& points) {
    LineLineIntersector intersector;
    vtkIdList* cellIds = vtkIdList::New();
    vtkIdList* ptIds = vtkIdList::New();

    std::vector<double> v0(3);
    std::vector<double> v1(3);

    // vector from start to finish
    double dp[] = {pEnd[0] - pBeg[0], pEnd[1] - pBeg[1], pEnd[2] - pBeg[2]};

    // find all the cells intersected by the line
    this->locator->FindCellsAlongLine((double*) &pBeg[0], 
                                      (double*) &pEnd[0], 
                                      this->tol, cellIds);

    //
    // collect the intersection points
    //

    // iterate over the cells along the line
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {

        // this cell Id
        vtkIdType cId = cellIds->GetId(i);

        // vertices, ptIds.GetNumberOfIds() should return 4
        // since we're dealing with quads only
        this->grid->GetCellPoints(cId, ptIds);


        // iterate over the quads' edges
        for (vtkIdType j0 = 0; j0 < 4; ++j0) {

            vtkIdType j1 = (j0 + 1) % 4;

            this->grid->GetPoint(ptIds->GetId(j0), &v0[0]);
            this->grid->GetPoint(ptIds->GetId(j1), &v1[0]);

            // look for an intersection
            intersector.setPoints(&pBeg[0], &pEnd[0], &v0[0], &v1[0]);

            if (! intersector.hasSolution(this->eps)) {
                // skip if no solution. FindCellsAlongLine may be too generous with
                // returning the list of intersected cells
                continue;
            }

            // we have a solution but it could be degenerate

            if (std::abs(intersector.getDet()) > this->eps) {
                // normal intersection, 1 solution
                std::vector<double> sol = intersector.getSolution();
                double lambRay = sol[0];
                double lambEdg = sol[1];

                // is it valid? Intersection must be within (p0, p1) and (q0, q1)
                if (lambRay >= (0. - this->eps100) && lambRay <= (1. + this->eps100)  && 
                    lambEdg >= (0. - this->eps100) && lambEdg <= (1. + this->eps100)) {

                    // compute the intersection point
                    double p[] = {pBeg[0] + lambRay*dp[0], pBeg[1] + lambRay*dp[1]};

                    // add to list
                    cIds.push_back(cId);
                    lambRays.push_back(lambRay);
                    points.push_back(std::vector<double>(p, p + 2)); // copies
                }
            }
            else {
                // det is almost zero
                // looks like the two lines (p0, p1) and (q0, q1) are overlapping
                // add the starting/ending points
                const std::pair<double, double> sol = intersector.getBegEndParamCoords();
                // linear param coord along line
                double lama = sol.first;
                double lamb = sol.second;
                // compute the points
                double pa[] = {pBeg[0] + lama*dp[0], pBeg[1] + lama*dp[1]};
                double pb[] = {pBeg[0] + lamb*dp[0], pBeg[1] + lamb*dp[1]};

                // add to lists both points
                cIds.push_back(cId);
                lambRays.push_back(lama);
                points.push_back(std::vector<double>(pa, pa + 2));

                cIds.push_back(cId);
                lambRays.push_back(lamb); // same Id as before
                points.push_back(std::vector<double>(pb, pb + 2));

            }

        } // end of edge loop

    } // end of cell loop

    // clean up
    cellIds->Delete();
    ptIds->Delete();
  
}


void 
PolysegmentIter::__collectLineGridSegments(const double p0[], const double p1[]) {

    // things we need to define
    vtkIdList* ptIds = vtkIdList::New();
    vtkGenericCell* cell = vtkGenericCell::New();
    vtkIdList* cellIds = vtkIdList::New();

    double xi[] = {0., 0., 0.};
    double closestPoint[] = {0., 0., 0.};
    double weights[8];

    int subId;
    double dist;
    std::vector<double> point(3, 0.0);
    
    // VTK wants 3d positions
    double pBeg[] = {p0[0], p0[1], 0.};
    double pEnd[] = {p1[0], p1[1], 0.};

    // add starting point
    vtkIdType cId = this->locator->FindCell(pBeg, this->eps, cell, xi, weights);
    if (cId >= 0) {
        // success
        this->cellIds.push_back(cId);
        this->xis.push_back(std::vector<double>(xi, xi + 2)); // copy
        this->ts.push_back(0.); // start of line
    }

    //
    // find all intersection points in between
    //

    std::vector<vtkIdType> cIds;
    std::vector<double> lambRays;
    std::vector< std::vector<double> > points;
    this->__collectIntersectionPoints(pBeg, pEnd, cIds, lambRays, points);

    // find the cell Id of the neighbouring cells
    size_t nXPts = cIds.size();
    for (size_t i = 0; i < nXPts; ++i) {

        vtkIdType cId = cIds[i];
        double lambRay = lambRays[i];
	// need to copy because points are 2-tuples and VTK always works with 3-tuples
        point.assign(&points[i][0], &points[i][2]);

	vtkCell* cell = this->grid->GetCell(cId);
        int found = cell->EvaluatePosition((double*) &point[0], closestPoint, 
                                                                subId, xi, dist, weights);
        if (found) {
            this->cellIds.push_back(cId);
            this->xis.push_back(std::vector<double>(xi, xi + 2));
            this->ts.push_back(lambRay);
        }
        else {
            std::cerr << "Warning: param coord search failed for point " << point[0] << ", " << point[1] 
                                                                         << " in cell " << cId << '\n';
        }
    }
 
    // add end point 
    cId = this->locator->FindCell(pEnd, this->eps, cell, xi, weights);
    if (cId >= 0) {
        // success
        this->cellIds.push_back(cId);
        this->xis.push_back(std::vector<double>(xi, xi + 2));
        this->ts.push_back(1.); // end of line
    }

    // clean up
    ptIds->Delete();
    cell->Delete();
    cellIds->Delete();

}

