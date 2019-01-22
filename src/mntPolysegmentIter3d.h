#include <MvVector.h>
#include <mntLineTriangleIntersector.h>
#include <vtkUnstructuredGrid.h>
//#include <vtkOBBTree.h>
#include <vtkCellLocator.h>
#include <map>
#include <algorithm>

#ifndef MNT_POLYSEGMENT_ITER_3D
#define MNT_POLYSEGMENT_ITER_3D

class PolysegmentIter3d {

public:

    /**
     * Constructor
     * @param grid instance of vtkUnstructuredGrid
     * @param locator vtkCellLocator instance attached to the above grid
     * @param p0 start point
     * @param p1 end point
     */
    PolysegmentIter3d(vtkUnstructuredGrid* grid, vtkCellLocator* locator, 
                      const double p0[], const double p1[]);

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
    const Vector<double>& getBegCellParamCoord() const;


    /**
     * Get the current cell parametric coordinates at the end of segment
     * @return 2d array
     */
    const Vector<double>& getEndCellParamCoord() const;


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

    // cell Ids for each intersection point
    std::vector<vtkIdType> cellIds;

    // begin/end parametric coordinate
    std::vector<double> segTas;
    std::vector<double> segTbs;

    // begin/end cell param coords
    std::vector< Vector<double> > segXias;
    std::vector< Vector<double> > segXibs;    

    vtkUnstructuredGrid* grid;

    //vtkOBBTree* locator;
    vtkCellLocator* locator;

    double eps;
    double eps100;
    size_t index;
};

#endif // MNT_POLYSEGMENT_ITER_3D
