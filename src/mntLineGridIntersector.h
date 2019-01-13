#include <vtkUnstructuredGrid.h>
#include <vtkOBBTree.h>
#include <MvVector.h>
#include <vector>

#ifndef MNT_LINE_GRID_INTERSECTOR
#define MNT_LINE_GRID_INTERSECTOR

class LineGridIntersector {

public:

    /**
     * Constructor
     * @param grid instance of vtkUnstructuredGrid
     */
    LineGridIntersector(vtkUnstructuredGrid* grid);

    /**
     * Destructor
     */
    ~LineGridIntersector();

    /**
     * Set start/end points
     */
    void setLine(const double pa[], const double pb[]);

    /**
     * Get the line parametric coordinates of intersections
     * @return value
     */
    const std::vector<double>& getIntersectionLineParamCoords() const;

    /**
     * Get the current cell parametric coordinates at the beginning of segment
     * @return 2d array
     */
    std::vector< Vector<double> > getIntersectionPoints() const;


private:

    // cell locator
    vtkOBBTree* locator;

    // start point
    Vector<double> pA;

    // direction of the line (ie pb - pa)
    Vector<double> direction;

    // parametric line coordinates of the intersections
    std::vector<double> tValues;

    // small tolerance 
    double tol;

};

#endif // MNT_LINE_GRID_INTERSECTOR
