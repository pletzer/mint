#include <limits> 
#include <vector>
#include <map>
#include <set>
#include <mntVecN.h>
#include <vtkVersion.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>

#ifndef VMT_CELL_LOCATOR
#define VMT_CELL_LOCATOR

/**
 * A class to quickly find a cell in an unstructured grid (as a drop-in alternative to vtkCellLocator)
 */

class vmtCellLocator {

public:

    static vmtCellLocator* New() {
        return new vmtCellLocator();
    }

    /**
     * Constructor
     */
    vmtCellLocator();

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
    void SetDataSet(vtkUnstructuredGrid* grid);

    /** 
     * Set average number of cells/faces per bucket
     * @param avgNumFacesPerBucket number
     */
    void SetNumberOfCellsPerBucket(int avgNumFacesPerBucket);

    /**
     * Build the locator
     */
    void BuildLocator();

    /**
     * Find cell given a target point
     * @param point target
     * @param tol2 tolerance
     * @param cell pointer to the cell
     * @param pcoords parametric coordinates of x in the cell (output)
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
     * Set the periodicity length in x
     * @param periodX length (0 means not periodic)
     * @note this can be set after calling BuildLocator
     */
    void setPeriodicityLengthX(double periodX);

    /**
     * Get the periodicity length in x
     * @return 0 if non-periodic, periodicity in first coordinate otherwise
     */
    double getPeriodicityLengthX() const;

    /**
     * Enable folding across poles
     */
    void enableFolding();

    /**
     * Declare whether the grid is a gnomonic (central-projection) cubed sphere, e.g. FV3's
     * @param isCubedSphere true if the grid is a cubed sphere
     * @note only meaningful together with a periodic-longitude-in-degrees grid (periodX
     *       set > 0); for such a grid, a straight (lon, lat) chord between two cell
     *       corners is an excellent approximation of that cube-face edge everywhere
     *       except close to a pole -- where, since gnomonic projection maps straight
     *       lines on a cube face to great circles, the actual edge is a great circle that
     *       a straight (lon, lat) chord can badly misrepresent (see
     *       invertSphericalBilinearPatch). A plain (possibly rotated) lon-lat grid has no
     *       such cube faces -- its cell edges (lines of constant longitude or latitude)
     *       are exactly what a straight (lon, lat) chord already represents (a meridian
     *       is itself a great circle; a circle of latitude is not, but is still not a
     *       cubed-sphere edge either) -- so leave this false (the default) for those.
     */
    void setCubedSphere(bool isCubedSphere) {
        this->isCubedSphere = isCubedSphere;
    }

    /**
     * Find all intersection points between line and the grid
     * @param pBeg start point of the line
     * @param pEnd end point of the line
     * @return list of (cellId, [lambda0, lambda1, periodXOffset]) pairs
     * @note lambda0/lambda1 are the linear parametric coordinates of the entry/exit points into/from the cell
     * @note periodXOffset is the periodic offset to add to pBeg[0] and pEnd[0]
     */
    std::vector< std::pair<vtkIdType, Vec4> >
    findIntersectionsWithLine(const Vec3& pBeg, const Vec3& pEnd);

    /**
     * Check if a point is inside a face
     * @param faceId face/cell Id
     * @param point point
     * @param tol tolerance
     * @return true if inside, false otherwise
     */
    bool containsPoint(vtkIdType faceId, const double point[3], double tol) const;

    /**
     * Check if a point is inside a face
     * @param faceId face/cell Id
     * @param point point
     * @param tol tolerance
     * @return true if inside, false otherwise
     * @note this version takes into account the the multiplicity of longitudes and the folding for points beyond 
     *       +-90 degrees
     */
    bool containsPointMultiValued(vtkIdType faceId, const double point[3], double tol) const;

    /**
     * Print the bucket to face indices map
     */
    void printBuckets() const;


protected:

private:

    vtkUnstructuredGrid* grid;

    // mid point of the longitudes
    double lambdaMid;

    // domain range
    double xmin[3];
    double xmax[3];

    double weights[8];

    std::vector<double> modPeriodX;

    // number of buckets in X and Y
    int numBucketsX;
    int numBucketsY;

    // maps a bucket to a list of faces
    std::map<int, std::set<vtkIdType> > bucket2Faces;

    // vector of one ({0}) or two components ({0, 1}). Use {0, 1} for
    // folding of over the poles. 
    std::vector<int> kFolding;

    // periodicity in x
    double periodX;

    // true if the grid is a gnomonic cubed sphere (see setCubedSphere)
    bool isCubedSphere;

    // for a cubed-sphere grid, every face's 4 corners as unit vectors in Cartesian XYZ
    // (flat array, faceId*4 + corner), computed once in BuildLocator instead of on every
    // containsPointCubedSphere call -- the same handful of candidate faces per bucket get
    // tested against many different target points, so this easily pays for itself
    std::vector<Vec3> cubedSphereVerts;

    // for a cubed-sphere grid, every face's centroid (mean of its 4 corners' unit
    // vectors, not itself a unit vector) and squared radius (largest squared distance
    // from that centroid to any of the 4 corners) -- the cheap pre-filter in
    // containsPointCubedSphere compares against these instead of recomputing them
    std::vector<Vec3> cubedSphereCentroid;
    std::vector<double> cubedSphereRadius2;

    /**
     * Adjust the longitude and latitude to account for the folding at the pole
     * @param point lon, lat in input and transformed lon, lat on output
     */
    inline void foldAtPole(double point[]) const {
        double sgnLambda = point[0] >= this->lambdaMid? 1.: -1.;
        double sgnTheta = point[1] >= 0? 1.: -1.;
        point[0] -= sgnLambda*180.;
        point[1] = sgnTheta*180 - point[1];
    }

    /**
     * Convert a (lon, lat) point, in degrees, to a unit vector in 3D Cartesian space
     * @param lonLatDeg pointer to lon, lat, in degrees
     * @return unit vector on the sphere
     */
    inline Vec3 lonLatDegToXYZ(const double lonLatDeg[2]) const {
        const double deg2rad = M_PI / 180.0;
        double lam = lonLatDeg[0] * deg2rad;
        double the = lonLatDeg[1] * deg2rad;
        double cosThe = std::cos(the);
        Vec3 xyz;
        xyz[0] = cosThe * std::cos(lam);
        xyz[1] = cosThe * std::sin(lam);
        xyz[2] = std::sin(the);
        return xyz;
    }

    /**
     * Spherical linear interpolation between two unit vectors
     * @param u start unit vector
     * @param w end unit vector
     * @param t interpolation parameter, in [0, 1]
     * @return unit vector, u at t=0, w at t=1, along the great-circle arc between them
     */
    inline Vec3 slerp(const Vec3& u, const Vec3& w, double t) const {
        double cosOmega = std::max(-1.0, std::min(1.0, dot(u, w)));
        double omega = std::acos(cosOmega);
        if (omega < 1.e-9) {
            // u and w (nearly) coincide -- avoid the 0/0 below, any t gives the same point
            return u;
        }
        double s = std::sin(omega);
        return (std::sin((1. - t)*omega)/s)*u + (std::sin(t*omega)/s)*w;
    }

    /**
     * Map parametric coordinates to a point on the exact spherical quad spanned by a
     * cubed-sphere face's 4 corners (a "spherical bilinear patch": interpolate along the
     * 0->1 and 3->2 edges at xsi, then between those two points at eta -- all by great
     * circle, via slerp)
     * @param xsi xsi parametric coordinate, in [0, 1] inside the face
     * @param eta eta parametric coordinate, in [0, 1] inside the face
     * @param verts the face's 4 corners, unit vectors in Cartesian XYZ, in the same
     *              0->1->2->3 corner order used everywhere else in this class
     * @return unit vector on the sphere
     */
    inline Vec3 sphericalBilinearMap(double xsi, double eta, const Vec3 verts[4]) const {
        Vec3 a = this->slerp(verts[0], verts[1], xsi);
        Vec3 b = this->slerp(verts[3], verts[2], xsi);
        Vec3 r = this->slerp(a, b, eta);
        return r / std::sqrt(dot(r, r));
    }

    /**
     * Find the parametric coordinates (xsi, eta) at which the spherical bilinear patch
     * spanned by a cubed-sphere face's 4 corners passes through a target point, by
     * Gauss-Newton iteration (finite-difference Jacobian)
     * @param target target point, unit vector in Cartesian XYZ
     * @param verts the face's 4 corners, unit vectors in Cartesian XYZ
     * @param xsi xsi parametric coordinate (output), whether or not the iteration converges
     * @param eta eta parametric coordinate (output), whether or not the iteration converges
     * @return true if the iteration converged
     * @note convergence to (xsi, eta) outside [0, 1] is expected and correct for a target
     *       point outside the face -- the values stay bounded and well behaved (verified
     *       against faces spanning up to a full pole-adjacent cubed-sphere panel corner),
     *       they are just not a valid location inside this particular face
     */
    bool invertSphericalBilinearPatch(const Vec3& target, const Vec3 verts[4],
                                       double& xsi, double& eta) const;

    /**
     * Check if a point is inside a cubed-sphere face, using the exact spherical bilinear
     * patch spanned by its 4 corners for both the containment decision and the
     * parametric coordinates -- see invertSphericalBilinearPatch
     * @param faceId face/cell Id
     * @param point point, (lon, lat) in degrees
     * @param tol tolerance on xsi, eta
     * @param pcoords parametric coordinates (output): pcoords[0] = xsi, pcoords[1] = eta,
     *                pcoords[2] = 0; valid whether or not the point turns out to be inside
     * @param weights bilinear interpolation weights from (xsi, eta) (output), one per
     *                corner, in the same 0->1->2->3 order
     * @return true if inside
     * @note deriving both the containment decision and the parametric coordinates from
     *       the same spherical model, instead of pairing an accurate containment test
     *       with vtkQuad's flat, straight-(lon,lat)-chord EvaluatePosition, is the point:
     *       the two can never disagree with each other close to a pole, where they used
     *       to -- a point the flat model sees as just outside its own corners can get
     *       assigned wildly extrapolated (xsi, eta), and from there arbitrarily large
     *       interpolation weights.
     */
    inline bool containsPointCubedSphere(vtkIdType faceId, const double point[3], double tol,
                                          double pcoords[3], double weights[8]) const {
        // corners, centroid and radius are all precomputed once, in BuildLocator -- see
        // cubedSphereVerts
        const Vec3* verts = &this->cubedSphereVerts[faceId*4];
        Vec3 target = this->lonLatDegToXYZ(point);

        // cheap pre-filter: the locator's buckets are sized for a handful-of-flops
        // containment test, so a bucket can hold many candidate cells for every one
        // that actually contains a given point -- reject the (typically large) majority
        // that are nowhere close before paying for a multi-iteration Newton solve.
        // Compare squared Cartesian distance to this cell's own centroid against its own
        // "radius" (the farthest corner from that centroid), inflated by a generous
        // safety factor: this can only ever over-accept (falling through to the exact
        // solve below), never wrongly reject a point that's actually inside -- every
        // point of the spherical bilinear patch stays within the spherical convex hull
        // of its 4 corners, which in turn stays within their centroid's own radius, so
        // safetyFactor > 1 is already conservative; the extra margin is just insurance.
        const double safetyFactor = 2.0;
        Vec3 dt = target - this->cubedSphereCentroid[faceId];
        if (dot(dt, dt) > safetyFactor*safetyFactor*this->cubedSphereRadius2[faceId]) {
            pcoords[0] = pcoords[1] = pcoords[2] = -1.0; // clearly outside, not solved
            weights[0] = weights[1] = weights[2] = weights[3] = 0.0;
            return false;
        }

        double xsi = 0.5, eta = 0.5;
        bool converged = this->invertSphericalBilinearPatch(target, verts, xsi, eta);

        pcoords[0] = xsi;
        pcoords[1] = eta;
        pcoords[2] = 0.0;

        weights[0] = (1. - xsi)*(1. - eta);
        weights[1] = xsi*(1. - eta);
        weights[2] = xsi*eta;
        weights[3] = (1. - xsi)*eta;

        return converged && xsi >= -tol && xsi <= 1. + tol && eta >= -tol && eta <= 1. + tol;
    }


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
        for (std::size_t i = 0; i < 2; ++i) {
            x[i] = (point[i] - this->xmin[i]) / (this->xmax[i] - this->xmin[i]); // must have some thickness!
        }

        // bucket coordinates
        int m = (int) std::floor(this->numBucketsX * x[0] + eps);
        int n = (int) std::floor(this->numBucketsY * x[1] + eps);

        // make sure the bucket coordinates fit in the domain
        m = std::max(0, std::min(this->numBucketsX - 1, m));
        n = std::max(0, std::min(this->numBucketsY - 1, n));

        // return flat array index
        return m * this->numBucketsY + n;
    }

    /**
     * Get the bucket index coordinates
     * @param bucketId flat bucket Id
     * @param m index (output)
     * @param n index (output)
     */
    inline void getBucketIndices(int bucketId, int* m, int* n) const {
        *m = bucketId / this->numBucketsY;
        *n = bucketId % this->numBucketsY;
    }

    /**
     * Collect the intersection points between line and cell
     * @param cellId cell Id
     * @param pBeg start point of the line
     * @param direction direction of the line (pEnd = pBeg + direction)
     * @return array of line parameter coordinates in increasing order
     * @note expect either 0 (no intersection) or 2 values (intersection) to be returned.
     *       start/end points qualify as intersection if they fall into the cell
     */
    std::vector<double> collectIntersectionPoints(vtkIdType cellId, 
                                                  const Vec3& pBeg,
                                                  const Vec3& direction);

    /**
     * Get the nodal points of the face
     * @param faceId face Id
     * @return list of points
     */
    inline std::vector<Vec3> getFacePoints(vtkIdType faceId) const {

#if( (VTK_MAJOR_VERSION == 8 && VTK_MINOR_VERSION == 90) || (VTK_MAJOR_VERSION >= 9) )
        // Paraview 5.8.5, may need to make this more general
        const vtkIdType* ptIds;
#else
        vtkIdType* ptIds;
#endif
        vtkIdType npts;
        this->grid->GetCellPoints(faceId, npts, ptIds);
        vtkPoints* points = this->grid->GetPoints();
        std::vector<Vec3> res(npts);
        for (vtkIdType i = 0; i < npts; ++i) {
            vtkIdType idx = ptIds[i];
            double* p = points->GetPoint(idx);
            res[i] = Vec3(p);
        }

        return res;
    }

};


#endif // VMT_CELL_LOCATOR
