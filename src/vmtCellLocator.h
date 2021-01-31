#include <vector>
#include <map>
#include <set>
#include <mntVecN.h>
#include <vtkVersion.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <limits>


#ifndef VMT_CELL_LOCATOR
#define VMT_CELL_LOCATOR

/**
 * A class to quickly find a cell in an unstructured grid  (as a drop-in alternative to vtkCellLocator)
 */

class vmtCellLocator {

public:

    /**
     * Constructor
     */
    vmtCellLocator();

    static vmtCellLocator* New() {
        return new vmtCellLocator();
    }

    /**
     * Destructor
     */
    ~vmtCellLocator() {
    }


    void Delete() {
        delete this;
    }

    /**
     * Set the grid
     * @param grid vtkUnstructuredGrid object
     */
    void SetDataSet(vtkUnstructuredGrid* grid);

    /** 
     * Set average number of cells/faces per bucket
     * @param avgNumFacesPerBucket number
     */
    void SetNumberOfCellsPerBucket(int avgNumFacesPerBucket);

    /**
     * Build the locator
     */
    void BuildLocator();

    /**
     * Find cell given a target point
     * @param point target
     * @param tol2 tolerance
     * @param cell pointer to the cell
     * @param pcoords parametric coodinates of x in the cell (output)
     * @param weights interpolation weights of the point
     * @return cell Id if found, < 0 otherwise
     */
    vtkIdType FindCell(const double point[3], double tol2, vtkGenericCell *cell, double pcoords[3], double *weights);

    /**
     * Find all the cells intersected by line
     * @param p0 start point
     * @param p1 end point
     * @param tol2 tolerance
     * @param cellIds list of cell Ids
     */
    void FindCellsAlongLine(const double p0[3], const double p1[3], double tol2, vtkIdList *cellIds);

    /**
     * Find intersection point
     * @param p0 start point
     * @param p1 end point
     * @param tol2 tolerance
     * @param t linear parametric coordinate of the intersection point (output)
     * @param x intersection point
     * @param pcoords cell parametric coordinates of the intersection point (output)
     * @param subId sub-cell Id (not used)
     * @param cellId cell Id 
     * @param cell pointer to the cell
     * @return 1 if found, 0 otherwise
     */
    int IntersectWithLine(const double p0[3], const double p1[3], double tol, double &t, double x[3], 
                          double pcoords[3], int &subId, vtkIdType &cellId, vtkGenericCell *cell);


    /**
     * Set the periodicity length in x
     * @param periodX length (0 means not periodic)
     * @note this can be set after calling BuildLocator
     */
    void setPeriodicityLengthX(double periodX);

    /**
     * Find all intersection points between line and the grid
     * @param pBeg start point of the line
     * @param pEnd end point of the line
     * @return list of (cellId, [lambda0, lambda1, periodXOffset]) pairs
     * @note lambda0/lambda1 are the linear parametric coordiates of the entry/exit points into/from the cell
     * @note periodXOffset is the periodic offset to add to pBeg[0] and pEnd[0]
     */
    std::vector< std::pair<vtkIdType, Vec3> >
    findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd);

    /**
     * Check if a point is indide a face
     * @param faceId face/cell Id
     * @param point point
     * @param tol tolerance
     * @return true if inside, false otherwise
     */
    bool containsPoint(vtkIdType faceId, const double point[3], double tol) const;

    /**
     * Print the bucket to face indices map
     */
    void printBuckets() const;


protected:

private:

    vtkUnstructuredGrid* grid;

    // domain range
    double xmin[3];
    double xmax[3];

    double weights[8];

    std::vector<double> modPeriodX;

    // number of buckets in X and Y
    int numBucketsX;
    int numBucketsY;

    // maps a bucket to a list of faces
    std::map<int, std::set<vtkIdType> > bucket2Faces;


    /**
     * Get the flat array index of a bucket containing a given point
     * @param point point
     * @return index
     * @note assumes there is thickness in the domain
     *       will return index ven if the point is outside the domain
     */
    inline int getBucketId(const double point[3]) const {

        // required to make sure std::floor that does not return the 
        // next integer below if we're close to an integer
        const double eps = 10 * std::numeric_limits<double>::epsilon();

        // normalize
        double x[3];
        for (size_t i = 0; i < 2; ++i) {
            x[i] = (point[i] - this->xmin[i]) / (this->xmax[i] - this->xmin[i]); // must have some thickness!
        }

        // bucket coordinates
        int m = (int) std::floor(this->numBucketsX * x[0] + eps);
        int n = (int) std::floor(this->numBucketsY * x[1] + eps);

        // make sure the bucket coordinates fit in the domain
        m = std::max(0, std::min(this->numBucketsX - 1, m));
        n = std::max(0, std::min(this->numBucketsY - 1, n));

        // return flat array index
        return m * this->numBucketsY + n;
    }

    /**
     * Get the bucket index coordinates
     * @param bucketId flat bucket Id
     * @param m index (output)
     * @param n index (output)
     */
    inline void getBucketIndices(int bucketId, int* m, int* n) const {
        *m = bucketId / this->numBucketsY;
        *n = bucketId % this->numBucketsY;
    }

    /**
     * Collect the intersection points between line and cell
     * @param cellId cell Id
     * @param pBeg start point of the line
     * @param direction direction of the line (pEnd = pBeg + direction)
     * @return array of line parameter coordinates in increasing order
     * @note expect either 0 (no intersection) or 2 values (intersection) to be returned.
     *       start/end points qualify as intersection if they fall into the cell
     */
    std::vector<double> collectIntersectionPoints(vtkIdType cellId, 
                                                  const Vec3& pBeg,
                                                  const Vec3& direction);

    /**
     * Get the nodal points of the face
     * @param faceId face Id
     * @return list of points
     */
    inline std::vector<Vec3> getFacePoints(vtkIdType faceId) const {

#if(VTK_MAJOR_VERSION >= 8 && VTK_MINOR_VERSION == 90)
        // Paraview 5.8.5, may need to make this more general
        const vtkIdType* ptIds;
#else
        vtkIdType* ptIds;
#endif
        vtkIdType npts;
        this->grid->GetCellPoints(faceId, npts, ptIds);
        vtkPoints* points = this->grid->GetPoints();
        std::vector<Vec3> res(npts);
        for (vtkIdType i = 0; i < npts; ++i) {
            vtkIdType idx = ptIds[i];
            double* p = points->GetPoint(idx);
            res[i] = Vec3(p);
        }

        return res;
    }

};


#endif // VMT_CELL_LOCATOR
