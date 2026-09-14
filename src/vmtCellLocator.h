#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <mntVecN.h>
#include <vector>
#include <utility>

#ifndef VMT_CELL_LOCATOR
#define VMT_CELL_LOCATOR

/**
 * Abstract interface for mint's cell locators: given a target point, find
 * which cell (if any) of an unstructured grid contains it, plus that
 * point's parametric coordinates/interpolation weights within the cell.
 *
 * Two concrete implementations:
 *  - vmtLonLatCellLocator: periodic-in-longitude, optionally pole-folding,
 *    for (lon, lat[, elev=0]) grids (this is the ORIGINAL vmtCellLocator,
 *    renamed -- see that class's own docstring for why (x,y)-only bucket
 *    binning is a valid simplification there).
 *  - vmtXYZCellLocator: for a genuinely 3D-embedded (x, y, z) surface mesh
 *    with no periodic seam -- delegates to a standard vtkStaticCellLocator,
 *    which indexes on all 3 coordinates and has no periodicity assumptions
 *    to get wrong (see test_cell_locator_xyz.py for the bug this fixes:
 *    vmtLonLatCellLocator's (x,y)-only bucket grid could altogether miss a
 *    bona fide interior point of a curved cell whose corners happen to
 *    project to a different bucket than the point itself).
 *
 * mnt_vectorinterp_buildLocator picks the right concrete class
 * automatically (from whether periodX > 0). To use a locator you built
 * and configured yourself instead -- e.g. sharing one locator across
 * several VectorInterp/PolylineIntegral objects, or picking the concrete
 * type explicitly rather than relying on that auto-selection -- construct
 * either subclass directly and hand it to setLocator.
 *
 * Covers what VectorInterp needs (plain point-in-cell queries via FindCell)
 * plus what PolysegmentIter needs (findIntersectionsWithLine, tracing a
 * straight line/edge across however many cells it crosses -- see that
 * method's own docstring for how the two concrete classes' answers differ:
 * vmtLonLatCellLocator solves this exactly in the shared flat (lon, lat)
 * plane, vmtXYZCellLocator approximates it for a genuinely curved surface
 * by treating the line as a straight 3D chord). PolylineIntegral itself
 * (as opposed to PolysegmentIter, which it uses) needs considerably more
 * beyond that (multi-valued containment near a pole, bucket introspection,
 * ...) and so uses vmtLonLatCellLocator directly rather than through this
 * interface -- that extra machinery is inherently lon-lat-specific and has
 * no vmtXYZCellLocator equivalent (yet).
 */
class vmtCellLocator {

public:

    virtual ~vmtCellLocator() {}

    virtual void Delete() = 0;

    /**
     * Set the grid
     * @param grid vtkUnstructuredGrid object
     */
    virtual void SetDataSet(vtkUnstructuredGrid* grid) = 0;

    /**
     * Set average number of cells/faces per bucket
     * @param avgNumFacesPerBucket number
     */
    virtual void SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) = 0;

    /**
     * Build the locator
     */
    virtual void BuildLocator() = 0;

    /**
     * Find cell given a target point
     * @param point target
     * @param tol2 tolerance
     * @param cell pointer to the cell
     * @param pcoords parametric coordinates of x in the cell (output)
     * @param weights interpolation weights of the point
     * @return cell Id if found, < 0 otherwise
     */
    virtual vtkIdType FindCell(const double point[3], double tol2, vtkGenericCell *cell,
                               double pcoords[3], double *weights) = 0;

    /**
     * Set the periodicity length in x
     * @param periodX length (0 means not periodic)
     * @note meaningless for vmtXYZCellLocator (no periodic seam); calling
     *       this with periodX > 0 on one is a usage error
     */
    virtual void setPeriodicityLengthX(double periodX) = 0;

    /**
     * Get the periodicity length in x
     * @return 0 if non-periodic, periodicity in first coordinate otherwise
     */
    virtual double getPeriodicityLengthX() const = 0;

    /**
     * Enable folding across poles
     * @note only meaningful for vmtLonLatCellLocator
     */
    virtual void enableFolding() = 0;

    /**
     * Declare whether the grid is a gnomonic (central-projection) cubed sphere
     * @param isCubedSphere true if the grid is a cubed sphere
     * @note only meaningful for vmtLonLatCellLocator
     */
    virtual void setCubedSphere(bool isCubedSphere) = 0;

    /**
     * Find all intersection points between a line and the grid -- used by
     * PolysegmentIter to break a (destination) edge into the sequence of
     * (source) cells it crosses.
     * @param pBeg start point of the line
     * @param pEnd end point of the line
     * @return list of (cellId, [lambda0, lambda1, periodXOffset, fold]) pairs
     * @note lambda0/lambda1 are the linear parametric coordinates of the
     *       entry/exit points into/from the cell, along pBeg->pEnd
     * @note periodXOffset/fold are always 0 for vmtXYZCellLocator (no
     *       periodicity or pole-folding for a plain embedded mesh); see
     *       vmtLonLatCellLocator.h for what they mean there
     */
    virtual std::vector< std::pair<vtkIdType, Vec4> >
    findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd) = 0;

};

#endif // VMT_CELL_LOCATOR
