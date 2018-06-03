#include <mntPolysegmentIter.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <MvVector.h>

PolysegmentIter::PolysegmentIter(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                                 const double p0[], const double p1[]) {

    // small tolerances 
    this->eps = 1.73654365e-14;
    this->eps100 = 100. * this->eps;
    this->tol = 1.e-3; // to determine if a point is inside a cell


    this->grid = grid;
    this->locator = locator;

    // cellIds, xis and ts are output
    this->cellIds.resize(0);
    this->xis.resize(0);
    this->ts.resize(0);
    this->__collectLineGridSegments(p0, p1);

    // re-arrange the data cellId -> [indx0, indx1, ...]
    // indx is index in the cellIds, xis and ts arrays
    std::map< vtkIdType, std::vector<size_t> > c2Inds;
    for (size_t i = 0; i < this->cellIds.size(); ++i) {
        vtkIdType cId = this->cellIds[i];
        std:map< vtkIdType, std::vector<size_t> >::iterator it = c2Inds.find(cId);
        if (it != c2Inds.end()) {
            // push_back
            it->second.push_back(i);
        }
        else {
            // create new entry
            std::vector<size_t> index(1, i);
            std::pair< vtkIdType, std::vector<size_t> > p(cId, index);
            c2Inds.insert(p);
        }
    }

    this->segCellIds.resize(0);
    this->segTas.resize(0);
    this->segTbs.resize(0);
    this->segXias.resize(0);
    this->segXibs.resize(0);
    this->segCoeffs.resize(0);
    // sort by t values
    TCmpFunctor f(this->ts);
    for (std::map< vtkIdType, std::vector<size_t> >::const_iterator it = c2Inds.begin();
        it != c2Inds.end(); ++it) {
        std::sort(it->second.begin(), it->second.end(), f);
        size_t n = it->second.size();
        for (size_t i = 0; i < n - 1; ++i) {
            size_t ia = it->second[i    ];
            size_t ib = it->second[i + 1];
            double ta = this->ts[ia];
            double tb = this->ts[ib];
            const std::vector<double>& xia = this->xis[ia];
            const std::vector<double>& xib = this->xis[ib];
            this->segCellIds.push_back(it->first);
            this->segTas.push_back(ta);
            this->segTbs.push_back(tb);
            this->segXias.push_back(xia);
            this->segXibs.push_back(xib);
            this->segCoeffs.push_back(1.0);
        }
    }

    size_t n = this->segCellIds.size();
    std::vector<size_t> inds(n);
    for (size_t i = 0; i < n; ++i)
        inds[i] = i;

    // sort arrays by ta values
    std::sort(inds.begin(), inds.end(), f);

    std::vector<vtkIdType> sCellIds = this->segCellIds;
    std::vector<double> sTas = this->segTas;
    std::vector<double> sTbs = this->segTbs;
    std::vector< std::vector<double> > sXias = this->segXias;
    std::vector< std::vector<double> > sXibs = this->segXibs;
    for (size_t i = 0; i < n; ++i) {
        size_t k = inds[i];
        this->segCellIds[i] = sCellIds[k];
        this->segTas[i] = sTas[k];
        this->segTbs[i] = sTbs[k];
        this->segXias[i] = sXias[k];
        this->segXibs[i] = sXibs[k];
    }

    // assign coefficients that account for duplicity, ie segments 
    // that are shared between two cells
    this->__assignCoefficientsToSegments();

    this->numSegs = this->segCellIds.size();

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
    this->index = -1;
}


void 
PolysegmentIter::next() {
    if (this->index < this->numSegs - 1) {
        this->index++;
    }
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
PolysegmentIter::getIndex() const {
    return this->index;
}


void
PolysegmentIter::__assignCoefficientsToSegments() {

    size_t n = this->segCellIds.size();

    // remove zero length segments
    for (int i = n - 1; i >= 0; --i) {
        double ta = this->segTas[i];
        double tb = this->segTbs[i];
        if (std::abs(tb - ta) < this->eps100) {
            this->segCellIds.erase(this->segCellIds.begin() + i);
            this->segTas.erase(this->segTas.begin() + i);
            this->segTbs.erase(this->segTbs.begin() + i);
            this->segXias.erase(this->segXias.begin() + i);
            this->segXibs.erase(this->segXibs.begin() + i);
            this->segCoeffs.erase(this->segCoeffs.begin() + i);
        }
    }

    // reduce contribution for overlapping segments. If two 
    // segments overlap then the coefficient of first segment
    // is set to 1.0 - overlap/(tb - ta). Assumes overlap 
    // can only happen for pairs of segment
    n = this->segCellIds.size();
    for (size_t i0 = 0; i0 < n - 1; ++i0) {
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

    double dp[] = {pEnd[0] - pBeg[0], pEnd[1] - pBeg[1], pEnd[2] - pBeg[2]};

    // find all the cells intersected by the line
    this->locator->FindCellsAlongLine((double*) &pBeg[0], 
                                      (double*) &pEnd[0], 
                                      this->tol, cellIds);

    // collect the intersection points in between
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {

        vtkIdType cId = cellIds->GetId(i);

        this->grid->GetCellPoints(cId, ptIds);

        // iterate over the quads' edges
        for (vtkIdType j0 = 0; j0 < 4; ++j0) {

            vtkIdType j1 = (j0 + 1) % 4;

            this->grid->GetPoint(ptIds->GetId(j0), &v0[0]);
            this->grid->GetPoint(ptIds->GetId(j1), &v1[0]);

            // look for an intersection
            intersector.setPoints(&pBeg[0], &pEnd[0], &v0[0], &v1[0]);
            if (! intersector.hasSolution(this->eps)) {
                continue;
            }

            if (std::abs(intersector.getDet()) > this->eps) {
                // normal intersection, 1 solution
                std::vector<double> sol = intersector.getSolution();
                double lambRay = sol[0];
                double lambEdg = sol[1];

                // is it valid? Intersection must be within (p0, p1) and (q0, q1)
                if (lambRay >= (0. - this->eps100) && lambRay <= (1. + this->eps100)  && 
                    lambEdg >= (0. - this->eps100) && lambEdg <= (1. + this->eps100)) {

                    double p[] = {pBeg[0] + lambRay*dp[0], pBeg[1] + lambRay*dp[1]};
                    cIds.push_back(cId);
                    lambRays.push_back(lambRay);
                    points.push_back(std::vector<double>(p, p + 2)); // is this a copy?
                }
            }
            else {
                // det is almost zero
                // looks like the two lines (p0, p1) and (q0, q1) are overlapping
                // add the starting/ending points
                const std::pair<double, double> sol = intersector.getBegEndParamCoords();
                double lama = sol.first;
                double lamb = sol.second;
                double pa[] = {pBeg[0] + lama*dp[0], pBeg[1] + lama*dp[1]};
                double pb[] = {pBeg[0] + lamb*dp[0], pBeg[1] + lamb*dp[1]};
                cIds.push_back(cId);
                lambRays.push_back(lama);
                points.push_back(std::vector<double>(pa, pa + 2));
                cIds.push_back(cId);
                lambRays.push_back(lamb);
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
    double weights[] = {0., 0., 0., 0.};

    int subId;
    double dist;
    
    // VTK wants 3d positions
    double pBeg[] = {p0[0], p0[1], 0.};
    double pEnd[] = {p1[0], p1[1], 0.};

    // add starting point
    vtkIdType cId = this->locator->FindCell(pBeg, this->eps, cell, xi, weights);
    if (cId >= 0) {
        this->cellIds.push_back(cId);
        this->xis.push_back(std::vector<double>(xi, xi + 2));
        this->ts.push_back(0.);
    }

    //
    // find all intersection points in between
    //

    std::vector<vtkIdType> cIds;
    std::vector<double> lambRays;
    std::vector< std::vector<double> > points;
    this->__collectIntersectionPoints(pBeg, pEnd, cIds, lambRays, points);

    // find the cell id of the neighbouring cells
    size_t nXPts = cIds.size();
    for (size_t i = 0; i < nXPts; ++i) {

        vtkIdType cId = cIds[i];
        double lambRay = lambRays[i];
        const std::vector<double>& point = points[i];

        int found = this->grid->GetCell(cId)->EvaluatePosition((double*) &point[0], closestPoint, 
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
        this->cellIds.push_back(cId);
        this->xis.push_back(std::vector<double>(xi, xi + 2));
        this->ts.push_back(1.);
    }

    // clean up
    ptIds->Delete();
    cell->Delete();
    cellIds->Delete();

}

