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

    /**
     * Trace the straight 3D chord pBeg->pEnd across however many cells it
     * crosses -- see the .cpp for why this is NOT a ray-cast (a chord
     * between two points on a curved surface dips strictly inside it
     * everywhere except at the endpoints, so a literal ray-vs-surface
     * intersection test would find zero intermediate cells) and what it
     * does instead. periodXOffset/fold in the returned Vec4 are always 0.
     */
    std::vector< std::pair<vtkIdType, Vec4> >
    findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd) override;

private:

    vtkStaticCellLocator* locator;
    vtkGenericCell* scratchCell;
    vtkUnstructuredGrid* grid;

    /**
     * Flat (non-spherical) bilinear map, corners in the usual 0->1->2->3
     * convention -- same shape as vmtLonLatCellLocator::sphericalBilinearMap
     * but without the slerp (this is a flat 3D patch, not one constrained
     * to a sphere).
     */
    inline Vec3 bilinearMap(double xi, double eta, const Vec3 verts[4]) const {
        return (1. - xi) * (1. - eta) * verts[0] + xi * (1. - eta) * verts[1]
             + xi * eta * verts[2] + (1. - xi) * eta * verts[3];
    }

    /**
     * Find the (xi, eta) at which the flat bilinear patch spanned by a
     * cell's 4 corners passes closest to target, by Gauss-Newton iteration
     * (finite-difference Jacobian) -- same recipe, tolerances and
     * iteration count as vmtLonLatCellLocator::invertSphericalBilinearPatch,
     * just for a flat (not great-circle-constrained) patch. xi, eta are
     * BOTH the initial guess (in) and the result (out): findIntersectionsWithLine
     * warm-starts each call from the previous one's converged (xi, eta),
     * which matters here -- an unconditional (0.5, 0.5) restart (what
     * vtkCell::EvaluatePosition uses internally) was tried first and found,
     * empirically, to fail to converge reliably during bisection, silently
     * producing wrong "outside" verdicts and hence wrong (too-short)
     * segments.
     * @return true if the iteration converged
     */
    bool invertBilinearPatch(const Vec3& target, const Vec3 verts[4],
                              double& xi, double& eta) const;

    /**
     * Is point p inside cell cellId, to within parametric tolerance tol on
     * its (xi, eta) coordinates? xi, eta are the warm-start guess (in) and
     * the converged solution (out) -- see invertBilinearPatch.
     * @note deliberately NOT going through FindCell/the bucket search --
     *       this checks one SPECIFIC, already-known cell directly, which
     *       is what the bisection in findIntersectionsWithLine needs
     */
    bool pointIsInCell(vtkIdType cellId, const double p[3], double tol,
                        double& xi, double& eta) const;

};

#endif // VMT_XYZ_CELL_LOCATOR
