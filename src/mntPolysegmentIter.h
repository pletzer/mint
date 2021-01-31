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
     * @param periodX length of the x periodic domain size (0 = non-periodic)
     */
    PolysegmentIter(vtkUnstructuredGrid* grid, vmtCellLocator* locator, 
                    const double p0[], const double p1[], 
                    double periodX=0.0);

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
     * @param v point (in/out)
     */
    void  __makePeriodic(Vec3& v);


    void __assignCoefficientsToSegments();

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
    double totalT;

    // either 360, 2*pi or 0 (if not periodic)
    double periodX;

    size_t index;
    size_t numSegs;

};

#endif // MNT_POLYSEGMENT_ITER
