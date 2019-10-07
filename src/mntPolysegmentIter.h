#include <mntVecN.h>
#include <mntLineLineIntersector.h>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <map>
#include <algorithm>

#ifndef MNT_POLYSEGMENT_ITER
#define MNT_POLYSEGMENT_ITER

class PolysegmentIter {

public:

    /**
     * Constructor
     * @param grid instance of vtkUnstructuredGrid
     * @param locator vtkCellLocator instance attached to the above grid
     * @param p0 start point
     * @param p1 end point
     * @param periodicityLength length of the x periodic domain size (0 = non-periodic)
     */
    PolysegmentIter(vtkUnstructuredGrid* grid, vmtCellLocator* locator, 
                    const double p0[], const double p1[], 
                    double periodicityLength=0.0);

    /**
     * Get the integrated linear parametric coordinates
     * @return value
     */
    double getIntegratedParamCoord() const;

    /**
     * Reset the counter
     */
    void reset();


    /**
     * Next segment
     * @return true if valid, false otherwise
     */
    bool next();

    /**
     * Get number of segments
     * @return number
     */
    size_t getNumberOfSegments() const;



    /**
     * Get the current cell Id
     * @return index
     */
    vtkIdType getCellId() const;


    /**
     * Get the current cell parametric coordinates at the beginning of segment
     * @return 2d array
     */
    const Vec3& getBegCellParamCoord() const;


    /**
     * Get the current cell parametric coordinates at the end of segment
     * @return 2d array
     */
    const Vec3& getEndCellParamCoord() const;


    /**
     * Get the current cell parametric coordinates at the beginning of segment
     * @return number
     */
    double getBegLineParamCoord() const;
        

    /**
     * Get the current cell parametric coordinates at the end of segment
     * @return number
     */
    double getEndLineParamCoord() const;

    /**
     * Get the coefficient accounting for duplicates
     * @return coefficient
     */
    double getCoefficient() const;


private:

    /**
     * Add/remove periodicity length to fit in domain
     * @param vBeg start point (in/out)
     * @param vEnd end point (in/out)
     */
    void  makePeriodic(Vec3& vBeg, Vec3& vEnd) {

        // fix start/end points if they fall outside the domain and the domain is periodic
        if (this->xPeriodicity > 0.) {
            double xmin = this->grid->GetBounds()[0];
            double xmax = this->grid->GetBounds()[1];
            if (vBeg[0] < xmin) {
                std::cerr << "Warning: adding periodicity length " << this->xPeriodicity << 
                         " to start point " << vBeg << "\n";
                vBeg[0] += this->xPeriodicity;
                std::cerr << "Now start point is " << vBeg << '\n';
            }
            else if (vBeg[0] > xmax) {
                std::cerr << "Warning: subtracting periodicity length " << this->xPeriodicity << 
                         " from start point " << vBeg << "\n";
                vBeg[0] -= this->xPeriodicity;
                std::cerr << "Now start point is " << vBeg << '\n';
            } 
            if (vEnd[0] < xmin) {
                std::cerr << "Warning: adding periodicity length " << this->xPeriodicity << 
                         " to end point " << vEnd << "\n";
                vEnd[0] += this->xPeriodicity;
                std::cerr << "Now end point is " << vEnd << '\n';
            }
            else if (vEnd[0] > xmax) {
                std::cerr << "Warning: subtracting periodicity length " << this->xPeriodicity << 
                         " from end point " << vEnd << "\n";
                vEnd[0] -= this->xPeriodicity;
                std::cerr << "Now end point is " << vEnd << '\n';
            }
        }
    }


    void __assignCoefficientsToSegments();

    /**
     * Collect all the intersection points
     * @param pBeg starting point
     * @param pEnd end point
     * @param cIds cell Ids (output)
     * @param lambRay line paramatric coordinates (output)
     * @param points intersection points (output)
     */
    void __collectIntersectionPoints(const double pBeg[], 
                                     const double pEnd[],
                                     std::vector<vtkIdType>& cIds,
                                     std::vector<double>& lambRays,
                                     std::vector<Vec2>& points);

    /**
     * Collect and store all the line-grid intersection points
     * @param p0 starting point of the line
     * @param p1 end point of the line 
     */
    void __collectLineGridSegments(const double p0[],
                                   const double p1[]);


    // cell Ids for each intersection point
    std::vector<vtkIdType> cellIds;
    
    // 2d cell parametric coordinates for each intersection point
    std::vector<Vec3> xis;

    // 1d line parametric coordinates for each intersection point
    std::vector<double> ts;

    // grid cell Ids for each segment
    std::vector<vtkIdType> segCellIds;

    // starten parametric line coordinates
    std::vector<double> segTas;
    std::vector<double> segTbs;

    // start/end cell parametric coordinates
    std::vector<Vec3> segXias;
    std::vector<Vec3> segXibs;

    // duplicity coefficient
    std::vector<double> segCoeffs;

    vtkUnstructuredGrid* grid;

    vmtCellLocator* locator;

    double eps;
    double eps100;
    double tol;
    double totalT;

    // either 360, 2*pi or 0 (if not periodic)
    double xPeriodicity;

    size_t index;
    size_t numSegs;

};

#endif // MNT_POLYSEGMENT_ITER
