#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <mntVecN.h>
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
    std::vector<Vec3> getIntersectionPoints() const;

    /**
     * Get direction of line
     * @return vector
     */
    const Vec3& getDirection() const;

    /**
     * Get start point
     * @return vector
     */
    const Vec3& getStartPoint() const;


private:

    // cell locator
    vtkCellLocator* locator;

    // start point
    Vec3 pA;

    // direction of the line (ie pb - pa)
    Vec3 direction;

    // parametric line coordinates of the intersections
    std::vector<double> tValues;

    // small tolerance 
    double tol;

};

#endif // MNT_LINE_GRID_INTERSECTOR
