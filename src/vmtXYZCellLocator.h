#include <vtkStaticCellLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <vmtCellLocator.h>

#ifndef VMT_XYZ_CELL_LOCATOR
#define VMT_XYZ_CELL_LOCATOR

/**
 * Cell locator for a genuinely 3D-embedded (x, y, z) surface mesh -- real
 * Cartesian coordinates, no lon/lat, no periodic seam -- implemented as a
 * thin wrapper around vtkStaticCellLocator, which indexes on all 3
 * coordinates (unlike vmtLonLatCellLocator's (x,y)-only bucket grid, see
 * that class's docstring and test_cell_locator_xyz.py for the bug that
 * causes).
 *
 * There is no periodicity or pole-folding here -- setPeriodicityLengthX
 * and setCubedSphere exist only to satisfy vmtCellLocator's shared
 * interface, and log an error if actually asked for (periodX > 0 /
 * isCubedSphere true): those only make sense for vmtLonLatCellLocator.
 */
class vmtXYZCellLocator : public vmtCellLocator {

public:

    static vmtXYZCellLocator* New() {
        return new vmtXYZCellLocator();
    }

    vmtXYZCellLocator();
    ~vmtXYZCellLocator();

    void Delete() override {
        delete this;
    }

    void SetDataSet(vtkUnstructuredGrid* grid) override;

    void SetNumberOfCellsPerBucket(int avgNumFacesPerBucket) override;

    void BuildLocator() override;

    /**
     * @note unlike vmtLonLatCellLocator, `cell` MAY be NULL here -- this
     *       class provides its own scratch vtkGenericCell in that case,
     *       since vtkStaticCellLocator::FindCell (unlike this class's own
     *       historical FindCell) requires a real one to write into.
     */
    vtkIdType FindCell(const double point[3], double tol2, vtkGenericCell *cell,
                       double pcoords[3], double *weights) override;

    void setPeriodicityLengthX(double periodX) override;

    double getPeriodicityLengthX() const override {
        return 0.0;
    }

    void enableFolding() override {
        // no poles to fold across for a plain embedded (x,y,z) mesh
    }

    void setCubedSphere(bool isCubedSphere) override;

private:

    vtkStaticCellLocator* locator;
    vtkGenericCell* scratchCell;

};

#endif // VMT_XYZ_CELL_LOCATOR
