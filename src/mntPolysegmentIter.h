#include <mntLineLineIntersector.h>
#include <vtkUnstructuredGrid.h>
#include <algorithm>

struct TCmpFunctor {
    TCmpFunctor(const std::vector<double>& ts) {
        this->tVals = ts;
    }
    bool operator()(size_t i, size_t j) {
        return (this->tVals[i] < this->tVals[j]);
    }
    std::vector<double> tVals;
}


struct PolysegmentIter {

    /**
     * Constructor
     * @param grid instance of vtkUnstructuredGrid
     * @param locator vtkCellLocator instance attached to the above grid
     * @param p0 start point
     * @param p1 end point
     */
    PolysegmentIter(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                    const double p0[], const double p1[], 
                    double t0=00, double t1=1.0) {

        // small tolerances 
        this->eps = 1.73654365e-14;
        this->eps100 = 100. * this->eps;
        this->tol = 1.e-3; // to determine if a point is inside a cell


        this->grid = grid;
        this->locator = locator;

        double dt = t1 - t0;

        // cellIds, xis and ts are output
        this->collectLineGridSegments(p0, p1);

        // re-arrange the data cellId -> [indx0, indx1, ...]
        // indx is index in the cellIds, xis and ts arrays
        std::map< vtkIdType, std::vector<size_t> > c2s;
        for (size_t i = 0; i < cellIds.size(); ++i) {
            vtkIdType cId = cellIds[i];
            const std::vector<double>& xi = this->xis[i];
            double t = this->ts[i];
            std:map< vtkIdType, std::vector<size_t> >::iterator it = c2s.find(cId);
            if (it != c2s.end()) {
                // append
                it->second.push_back(i);
            }
            else {
                // create new entry
                std::vector<size_t> v(1, i);
                std::pair< vtkIdType, std::vector<size_t> > p(cId, v);
                c2s.insert(p);
            }
        }


        // create list of [(i0, i1), ...] tuples
        this->segIndices.erase();
        for (std::map<vtkIdType, std::vector<size_t> >::const_iterator 
        	   it = c2s.begin(); it != c2s.end(); ++it) {
        	if (it->second.size() == 2) {
                // two points in each cell
                this->segIndices.push_back(it->second);
        	}
        }
        // sort by t values
        TCmpFunctor f(ts);
        std::sort(this->segIndices.begin(), this->segIndices.end(), f);

        // assign coefficients that account for duplicity, ie segments 
        // that are shared between two cells
        this->assignCoefficientsToSegments();

        this->numSegs = len(this->coeffs);

        this->totalT = 0;
        for (size_t i = 0; i < segIndices.size(); ++i) {
            size_t ia = segIndices[i].first;
            size_t ib = segIndices[i].second;
            double ta = ts[ia];
            double tb = ts[ib];
            this->totalT += (tb - ta) * coeffs[i];
        }            

        // reset the iterator
        this->reset()
    }


    /**
     * Get the integrated linear parametric coordinates
     * @return value
     */
    double getIntegratedParamCoord() const {
        return this->totalT;
    }


    /**
     * Reset the counter
     */
    void reset() {
        this->index = -1;
    }


    /**
     * Next segment
     */
    void next() {
        if (this->index < this->numSegs - 1) {
            this->index++;
        }
    }


    /**
     * Get the current cell Id
     * @return index
     */
    vtkIdType getCellId() const {

        return this->cellIds[this->index];
    }


    def getBegCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the beginning of segment
        @return 2d array
        """
        return this->segment[1][
        

    def getEndCellParamCoord(self):
        """
        Get the current cell parametric coordinates at the end of segment
        @return 2d array
        """
        return this->segment[1][2]
 

    def getBegLineParamCoord(self):
        """
        Get the current line parametric coordinates at the beginning of segment
        @return 2d array
        """
        return this->segment[0][0]
        

    def getEndLineParamCoord(self):
        """
        Get the current line parametric coordinates at the end of segment
        @return 2d array
        """
        return this->segment[0][1]

    def getCoefficient(self):
        """
        Get the coefficient accounting for duplicates
        @return coefficient
        """
        return this->segment[1][3]
 

    def getIndex(self):
        """
        Get the current index
        @return index
        """
        return this->index


    def __assignCoefficientsToSegments(self):

        n = len(this->data)

        // remove zero length segments
        for i in range(n - 1, -1, -1):
            s = this->data[i]
            ta, tb = s[0]
            if abs(tb - ta) < this->eps100:
                del this->data[i]

        // reduce contribution for overlapping segmmnts. If two 
        // segments overlap then the coefficient of first segment
        // is set to 1.0 - overlap/(tb - ta). Assumes overlap 
        // can only happen for pairs of segment
        n = len(this->data)
        for i0 in range(n - 1):
            i1 = i0 + 1
            s0 = this->data[i0]
            s1 = this->data[i1]
            ta0, tb0 = s0[0]
            ta1, tb1 = s1[0]
            overlap = max(0., min(tb0, tb1) - max(ta1, ta0))
            s0[1][3] = 1.0 - overlap/(tb0 - ta0)


    def collectIntersectionPoints(self, pBeg, pEnd):
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


    def __collectLineGridSegments(self, p0, p1):
        """
        Collect all the line-grid intersection points
        @param p0 starting point of the line
        @param p1 end point of the line 
        @return list of [ (cellId, xi, t), ...]
        """

        res = []

        // things we need to define
        ptIds = vtk.vtkIdList()
        cell = vtk.vtkGenericCell()
        cellIds = vtk.vtkIdList()
        subId = vtk.mutable(-1)
        dist = vtk.mutable(0.)
        xi = numpy.zeros((3,), numpy.float64)
        pBeg = numpy.zeros((3,), numpy.float64)
        pEnd = numpy.zeros((3,), numpy.float64)
        point = numpy.zeros((3,), numpy.float64)
        closestPoint = numpy.zeros((3,), numpy.float64)
        weights = numpy.array((4,), numpy.float64)

        // VTK wants 3d positions
        pBeg[:] = p0[0], p0[1], 0.0
        pEnd[:] = p1[0], p1[1], 0.0

        deltaPos = pEnd - pBeg

        // add starting point
        cId = this->locator.FindCell(pBeg, this->eps, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 0.) )
        else:
            pass
            #print('Warning: starting point {} not found!'.format(p0))

        tLast = 0.0
        cIdLast = cId
        fullSegment = False

        #
        // find all intersection points in between
        #

        // let's be generous with the collection of cells
        intersections = this->__collectIntersectionPoints(pBeg, pEnd)

        // find the cell id of the neighbouring cells
        for cId, lambRay, point in intersections:

            found = this->grid.GetCell(cId).EvaluatePosition(point, closestPoint, subId, xi, dist, weights)
            if found:
                res.append( (cId, xi[:2].copy(), lambRay) )
            else:
                print('Warning: param coord search failed point {} in cell {}'.format(point, cId))

            
        // add last point 
        cId = this->locator.FindCell(pEnd, this->eps, cell, xi, weights)
        if cId >= 0:
            res.append( (cId, xi[:2].copy(), 1.) )
        else:
            pass
            #print('Warning: end point {} not found!'.format(p1))

        return res

    // cell Ids for each intersection point
    std::vector<vtkIdType> cellIds;
    
    // 2d cell parametric coordinates for each intersection point
    std::vector< std::vector<double> > xis;

    // 1d line parametric coordinates for each intersection point
    std::vector<double> ts;

    // array of coefficients, one for each segment
    std::vector<double> coeffs;

    // create list of [(i0, i1), ...] tuples
    std::vector< std::pair<size_t, size_t> > segIndices;


};
