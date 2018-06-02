#include <PolysegmentIter.h>

PolysegmentIter::PolysegmentIter(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                                 const double p0[], const double p1[],) {

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
            // append
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
        this->segCellIds[i] = sCellids[k];
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
PolysegmentIter::getCoefficient() {
    return this->segCoeffs[this->index];
}
 

size_t
PolysegmentIter::getIndex() {
    return this->index;
}


void
PolysegmentIter::__assignCoefficientsToSegments() {

    size_t n = this->segCellIds.size();

    // remove zero length segments
    for (size_t i = n - 1; i >= -1, --i) {
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
                                             const double pEnd[]) {
        """
        Collect all the intersection points
        @param pBeg starting point
        @param pEnd end point
        @return [(cellId, lambda, point), ...]
        @note lambda is the linear parametric coordinate along the line
        """

        res = []

        intersector = LineLineIntersector()
        cellIds = vtk.vtkIdList()
        ptIds = vtk.vtkIdList()

        dp = pEnd - pBeg

        // find all the cells intersected by the line
        this->locator.FindCellsAlongLine(pBeg, pEnd, this->tol, cellIds)

        // collect the intersection points in between
        for i in range(cellIds.GetNumberOfIds()):

            cId = cellIds.GetId(i)

            this->grid.GetCellPoints(cId, ptIds)

            // iterate over the quads' edges
            for j0 in range(4):

                j1 = (j0 + 1) % 4
                v0 = numpy.array(this->grid.GetPoint(ptIds.GetId(j0)))
                v1 = numpy.array(this->grid.GetPoint(ptIds.GetId(j1)))

                // look for an intersection
                intersector.setPoints(pBeg[:2], pEnd[:2], v0[:2], v1[:2])
                if not intersector.hasSolution(this->eps):
                    continue

                if abs(intersector.getDet()) > this->eps:
                    // normal intersection, 1 solution
                    lambRay, lambEdg = intersector.getSolution()

                    // is it valid? Intersection must be within (p0, p1) and (q0, q1)
                    if lambRay >= 0. - this->eps100 and lambRay <= 1. + this->eps100 and \
                        lambEdg >= 0. - this->eps100 and lambEdg <= 1. + this->eps100:

                        point = pBeg + lambRay*dp
                        res.append( (cId, lambRay, point) )

                else:
                    // det is almost zero
                    // looks like the two lines (p0, p1) and (q0, q1) are overlapping
                    // add the starting/ending points 
                    lama, lamb = intersector.getBegEndParamCoords()
                    pa = pBeg + lama*dp
                    pb = pBeg + lamb*dp
                    res.append( (cId, lama, pa) )
                    res.append( (cId, lamb, pb) )

        return res
}


void 
PolysegmentIter::__collectLineGridSegments(const double p0[],
                                           const double p1[]) {
    // things we need to define
    vtkIdList* ptIds = vtkIdList::New();
    vtkGenericCell* cell = vtkGenericCell::New();
    vtkIdList* cellIds = vtkIdList::New();
    double xi[] = {0., 0., 0.};
    double point[] = {0., 0., 0.};
    double closestPoint[] = {0., 0., 0.};
    double weights[] = {0., 0., 0., 0.};

    vtkIdType subId;
    double dist;
    
    // VTK wants 3d positions
    double pBeg[] = {p0[0], p0[1], 0.};
    double pEnd[] = {p1[0], p1[1], 0.};

    // add starting point
    vtkTypeId cId = this->locator.FindCell(pBeg, this->eps, cell, xi, weights);
    if (cId >= 0) {
        this->cellIds.push_back(cId);
        this->xis.push_back(std::vector<double>(xi, xi + 2));
        this->ts.push_back(0.);
    }

    //
    // find all intersection points in between
    //

    // let's be generous with the collection of cells
    intersections = this->__collectIntersectionPoints(pBeg, pEnd);

        # find the cell id of the neighbouring cells
        for cId, lambRay, point in intersections:

            found = this->grid.GetCell(cId).EvaluatePosition(point, closestPoint, subId, xi, dist, weights)
            if found:
                this->cellIds.append(cId)
                this->xis.append(xi[:2].copy())
                this->ts.append(lambRay)
            else:
                print('Warning: param coord search failed point {} in cell {}'.format(point, cId))

            
        # add last point 
        cId = this->locator.FindCell(pEnd, this->eps, cell, xi, weights)
        if cId >= 0:
            this->cellIds.append(cId)
            this->xis.append(xi[:2].copy())
            this->ts.append(1.)
        else:
            pass
            #print('Warning: end point {} not found!'.format(p1))

}

