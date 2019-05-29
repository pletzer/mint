#include <vector>
#include <map>
#include <set>
#include "MvVector.h"
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <limits>


#ifndef VMT_CELL_LOCATOR
#define VMT_CELL_LOCATOR

/**
 * A class to quickly find a cell in an unstructured grid  (as a drop-in alternative to vtlCellLocator)
 */

class vmtCellLocator {

public:

    /**
     * Constructor
     */
    vmtCellLocator() {
        this->grid = NULL;
        this->points = NULL;
        this->numBucketsX = 10;
        double big = std::numeric_limits<double>::max();
        for (size_t i = 0; i < 3; ++i) {
            this->xmin[i] = big;
            this->xmax[i] = -big;
        }
    }

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
    void SetDataSet(vtkUnstructuredGrid* grid) {
        this->grid = grid;
        this->points = grid->GetPoints();
        double* bounds = grid->GetBounds();
        this->xmin[0] = bounds[0];
        this->xmax[0] = bounds[1];
        this->xmin[1] = bounds[2];
        this->xmax[1] = bounds[3];
        this->xmin[2] = bounds[4];
        this->xmax[2] = bounds[5];
    }

    /** 
     * Set average number of cells/faces per bucket
     * @param avgNumFacesPerBucket number
     */
    void SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) {

        vtkIdType numFaces = this->grid->GetNumberOfCells();

        // number of buckets along one dimension (2D)
        this->numBucketsX = (int) std::max(1.0, 
                              std::sqrt((double) numFaces / (double) avgNumFacesPerBucket)
                                      );

    }

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
     * Find all intersection points between line and the grid
     * @param pBeg start point of the line
     * @param pEnd end point of the line
     * @return list of (cellId, [lambda0, lambda1]) pairs
     */
    std::vector< std::pair<vtkIdType, std::vector<double> > >
    findIntersectionsWithLine(const Vector<double>& pBeg, const Vector<double>& pEnd);

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
    void printBuckets() const {
        for (const auto& b2f : this->bucket2Faces) {
            int bucketId = b2f.first;
            int m, n;
            this->getBucketIndices(bucketId, &m, &n);
            std::cout << "bucket " << bucketId << " (" << m << ',' << n << ") contains faces ";
            for (const auto& faceId : b2f.second) {
	        std::cout << faceId << ' ';
            }
            std::cout << '\n';
        }
    }

protected:

private:

    vtkUnstructuredGrid* grid;

    vtkPoints* points;

    // domain range
    double xmin[3];
    double xmax[3];

    double weights[8];

    // number of buckets in X and Y
    int numBucketsX;

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
        int n = (int) std::floor(this->numBucketsX * x[1] + eps);

        // make sure the bucket coordinates fit in the domain
        m = std::max(0, std::min(this->numBucketsX - 1, m));
        n = std::max(0, std::min(this->numBucketsX - 1, n));

        // return flat array index
        return m * this->numBucketsX + n;
    }

    /**
     * Get the bucket index coordinates
     * @param bucketId flat bucket Id
     * @param m index (output)
     * @param n index (output)
     */
    inline void getBucketIndices(int bucketId, int* m, int* n) const {
        *m = bucketId / this->numBucketsX;
        *n = bucketId % this->numBucketsX;
    }

    /**
     * Collect the intersection points between line and cell
     * @param cellId cell Id
     * @param pBeg start point of the line
     * @param pEnd end point of the line
     * @return array of line parameter coordinates in increasing order
     * @note expect either 0 (no intersection) or 2 values (intersection) to be returned.
     *       start/end points qualify as intersection if they fall into the cell
     */
    std::vector<double> collectIntersectionPoints(vtkIdType cellId, 
                                                  const double pBeg[3],
                                                  const double pEnd[3]);

    /**
     * Get the nodal points of the face
     * @param faceId face Id
     * @return list of points
     */
    inline std::vector< Vector<double> > getFacePoints(vtkIdType faceId) const {

        vtkIdType* ptIds;
        vtkIdType npts;
        this->grid->GetCellPoints(faceId, npts, ptIds);

        std::vector< Vector<double> > res(npts);
        for (size_t i = 0; i < npts; ++i) {
            vtkIdType idx = ptIds[i];
            double* p = this->points->GetPoint(idx);
            res[i] = Vector<double>{p[0], p[1], p[2]};
        }

        return res;
    }

};


#endif // VMT_CELL_LOCATOR
