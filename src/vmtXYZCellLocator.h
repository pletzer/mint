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

    /**
     * Parametric coordinates of p in cell cellId by central projection, the
     * same map findIntersectionsWithLine uses (see invertBilinearPatch)
     */
    void getParametricCoords(vtkUnstructuredGrid* grid, vtkIdType cellId,
                             const Vec3& p, Vec3& pcoords) override;

private:

    vtkStaticCellLocator* locator;
    vtkGenericCell* scratchCell;
    vtkUnstructuredGrid* grid;

    // centre of the central projection of points onto the cells, see setProjectionCentre
    Vec3 projCentre;

    /**
     * Choose the projection centre: the origin if every cell faces it from
     * the same side (e.g. a sphere), otherwise a point far off the mesh
     * along its mean normal (e.g. a flat mesh). Each ray from the centre
     * must meet the mesh at most once.
     */
    void setProjectionCentre();

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
     * Find the (xi, eta) at which the ray from the projection centre through
     * target meets the flat bilinear patch spanned by a cell's 4 corners
     * (central projection), by Newton iteration on (xi, eta, lambda) with
     * bilinearMap(xi, eta) = lambda*target, target and verts relative to the
     * projection centre. Two cells sharing an edge then map a point whose
     * ray crosses that edge to the same point of the edge, so the cells tile
     * the surface with no gaps or overlaps, and line integrals telescope
     * across cell boundaries. A point on the patch itself (e.g. on a flat
     * mesh) maps to itself. Falls back to closestPointOnBilinearPatch if the
     * ray is (nearly) parallel to the patch. xi, eta are both the initial
     * guess (in) and the result (out).
     * @return false if the iteration did not converge or the patch is behind
     *         the centre (lambda <= 0, e.g. the antipodal cell of a sphere),
     *         true otherwise (including the fallback)
     */
    bool invertBilinearPatch(const Vec3& target, const Vec3 verts[4],
                              double& xi, double& eta) const;

    /**
     * Find the (xi, eta) at which the flat bilinear patch passes closest to
     * target, by Gauss-Newton iteration (finite-difference Jacobian). Not
     * consistent between neighbouring cells for points off the patches;
     * only used as a fallback by invertBilinearPatch.
     * @return true if the iteration converged
     */
    bool closestPointOnBilinearPatch(const Vec3& target, const Vec3 verts[4],
                                     double& xi, double& eta) const;

    /**
     * Get the 4 corners of a cell, relative to the projection centre
     */
    void getCellVerts(vtkIdType cellId, Vec3 verts[4]) const;

    /**
     * Unit normals of the 4 planes through the projection centre and the edges of a
     * cell, oriented towards the cell's interior (zero for a degenerate edge)
     */
    void getEdgePlanes(vtkIdType cellId, Vec3 normals[4]) const;

    /**
     * Is p inside the cone of a cell (on the inner side of its 4 edge
     * planes, to within relative tolerance eps)?
     */
    bool coneContains(const Vec3 normals[4], const Vec3& p, double eps) const;

    /**
     * Line parameter >= t0 at which pBeg + t*direction leaves a cell's
     * cone, capped at 1
     */
    double coneExit(const Vec3 normals[4], const Vec3& pBeg, const Vec3& direction,
                    double t0) const;

    /**
     * Find the cell that the line pBeg + t*direction is in just after t = t0,
     * with pBeg relative to the projection centre
     * @param tExit line parameter at which the line leaves that cell (output)
     * @return cell Id, -1 if there is none
     */
    vtkIdType findCellAlongLine(const Vec3& pBeg, const Vec3& direction, double t0,
                                double& tExit);

};

#endif // VMT_XYZ_CELL_LOCATOR
