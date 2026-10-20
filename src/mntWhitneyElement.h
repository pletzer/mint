#include "mntLIBRARY_API.h"
#include <limits> // required by vtkUnstructuredGrid
#include <vtkType.h>
#include <mntVecN.h>
#include <mntGrid.h>

#ifndef MNT_WHITNEY_ELEMENT
#define MNT_WHITNEY_ELEMENT

/** upper bound on the number of degrees of freedom per cell of any element */
#define MNT_MAX_DOFS_PER_CELL 12

/**
 * Abstract interface for an edge (1-form, H(curl)) element on the reference
 * quad (xi, eta) in [0, 1]^2.
 *
 * A field is represented in each cell as sum_i dofs[i] w_i, where the basis
 * functions w_i = a_i(xi, eta) dxi + b_i(xi, eta) deta are dual to the
 * degrees of freedom (line integrals over edges or sub-edges). (xi, eta) are
 * always the parametric coordinates of the grid cell cellId found by the
 * locator, even if the element spans several grid cells. The physical
 * vector follows from the covariant components (a, b) through
 *
 *    edge (W1) vector:  a * grad xi + b * grad eta
 *    face (W2) vector:  (b * dr/dxi - a * dr/deta) / jac
 *
 * Concrete implementations:
 *  - WhitneyElementOrder1: lowest order Whitney element, one degree of
 *    freedom per edge (the edge integral). This is what mint has always
 *    used.
 *  - WhitneyElementMacro2: second order element on 2x2 blocks of grid
 *    cells, whose degrees of freedom are the 12 edge integrals of the block
 *    (see mntWhitneyElementMacro2.h).
 *
 * VectorInterp, PolylineIntegral and RegridEdges evaluate and integrate
 * the basis exclusively through this interface, so that a higher order
 * element can be plugged in.
 */
class WhitneyElement {

public:

    virtual ~WhitneyElement() {}

    /**
     * @return polynomial order of the element (1 = lowest order)
     */
    virtual int getOrder() const = 0;

    /**
     * @return number of degrees of freedom per cell
     */
    virtual int getNumDofs() const = 0;

    /**
     * Evaluate the covariant reference components of each basis function
     * @param cellId grid cell Id
     * @param xi first parametric coordinate
     * @param eta second parametric coordinate
     * @param a dxi components, array of size getNumDofs() (output)
     * @param b deta components, array of size getNumDofs() (output)
     */
    virtual void evalBasis(vtkIdType cellId, double xi, double eta, double a[], double b[]) const = 0;

    /**
     * Compute the covariant reference components of a field
     * @param cellId grid cell Id
     * @param xi first parametric coordinate
     * @param eta second parametric coordinate
     * @param dofs degrees of freedom of the cell, array of size getNumDofs()
     * @param a dxi component (output)
     * @param b deta component (output)
     */
    virtual void getCovariantComponents(vtkIdType cellId, double xi, double eta, const double dofs[],
                                        double* a, double* b) const;

    /**
     * Integrate each basis function along a straight segment in parameter space
     * @param cellId grid cell Id, the segment must lie within this cell
     * @param xia start point of the segment, in parametric coordinates
     * @param xib end point of the segment, in parametric coordinates
     * @param w integrals, array of size getNumDofs() (output)
     */
    virtual void integrateBasis(vtkIdType cellId, const Vec3& xia, const Vec3& xib, double w[]) const = 0;

    /**
     * Get the sign that converts a degree of freedom from the default
     * orientation (edges pointing in the positive xi/eta direction) to the
     * counterclockwise orientation of the cell boundary
     * @param dofIndex local degree of freedom index
     * @return +1 or -1
     */
    virtual int getCounterclockwiseSign(int dofIndex) const = 0;

    /**
     * Map a degree of freedom to its location in the user's field data
     *
     * The field data always hold edge integrals of the grid, either cell by
     * cell (numCells*MNT_NUM_EDGES_PER_QUAD) or on unique edges (numEdges),
     * whatever the element. Weights are stored per (cellId, dofIndex), where
     * cellId is the grid cell returned by the locator; this method tells how
     * to fetch the value of that degree of freedom.
     *
     * @param grid grid
     * @param cellId cell Id
     * @param dofIndex local degree of freedom index, 0 <= dofIndex < getNumDofs()
     * @param placement MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA
     * @param index index into the field data (output)
     * @param sign orientation sign to apply to the data value (output)
     * @return error code (0 = OK)
     */
    virtual int getDofDataIndex(Grid_t* grid, vtkIdType cellId, int dofIndex, int placement,
                                std::size_t* index, int* sign) const = 0;

    /**
     * Gather the degrees of freedom of a cell, with orientation signs applied
     * @param grid grid
     * @param cellId cell Id
     * @param data field data, either cell by cell or on unique edges
     * @param placement MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA
     * @param dofs degrees of freedom, array of size getNumDofs() (output)
     * @return number of degrees of freedom that could not be found (0 = OK),
     *         these are set to zero
     */
    virtual int getCellDofs(Grid_t* grid, vtkIdType cellId, const double data[], int placement,
                            double dofs[]) const;
};


/**
 * Lowest order Whitney edge element on the quad, one degree of freedom per
 * edge. Edges are numbered 0 (south), 1 (east), 2 (north) and 3 (west) and
 * point in the positive xi or eta direction:
 *
 *    w0 = (1 - eta) dxi,  w1 = xi deta,  w2 = eta dxi,  w3 = (1 - xi) deta
 */
class WhitneyElementOrder1 : public WhitneyElement {

public:

    int getOrder() const override { return 1; }

    int getNumDofs() const override { return 4; }

    void evalBasis(vtkIdType cellId, double xi, double eta, double a[], double b[]) const override;

    void getCovariantComponents(vtkIdType cellId, double xi, double eta, const double dofs[],
                                double* a, double* b) const override;

    void integrateBasis(vtkIdType cellId, const Vec3& xia, const Vec3& xib, double w[]) const override;

    int getCounterclockwiseSign(int dofIndex) const override;

    int getDofDataIndex(Grid_t* grid, vtkIdType cellId, int dofIndex, int placement,
                        std::size_t* index, int* sign) const override;
};

/**
 * Create an edge element of a given order
 * @param order 1 for WhitneyElementOrder1, 2 for WhitneyElementMacro2
 * @param grid grid, needed for order 2 (the grid is tiled into 2x2 blocks of cells)
 * @param element the new element, owned by the caller (output). Left unchanged on error
 * @param biquadratic order 2 only: map each 2x2 block with its biquadratic 9-node map
 *                    (true) or use the piecewise bilinear cell coordinates (false)
 * @return error code (0 = OK)
 */
int mnt_whitneyelement_create(int order, Grid_t* grid, WhitneyElement** element, bool biquadratic = true);

/**
 * Replace an element by a new element of a given order
 * @param element the element, deleted and replaced on success, unchanged on error
 * @param order 1 or 2
 * @param grid grid, needed for order 2
 * @return error code (0 = OK)
 */
int mnt_whitneyelement_setOrder(WhitneyElement** element, int order, Grid_t* grid);

/**
 * Rebuild an element for a new grid, keeping its order. If the new grid
 * cannot hold an element of that order, the element is replaced by the
 * lowest order element
 * @param element the element
 * @param grid the new grid
 * @return error code (0 = OK, otherwise the element fell back to order 1)
 */
int mnt_whitneyelement_rebuild(WhitneyElement** element, Grid_t* grid);

/*
 * C API, to evaluate an edge element directly (e.g. from Python), outside of
 * VectorInterp, PolylineIntegral and RegridEdges
 */

/**
 * Constructor
 * @param self the new element (output)
 * @param order 1 (lowest order Whitney element) or 2 (second order element on 2x2
 *              blocks of grid cells)
 * @param grid grid, needed for order 2. The element borrows it, keep it alive
 * @param biquadratic order 2 only: 1 to map each 2x2 block with its biquadratic
 *                    9-node map, 0 to use the piecewise bilinear cell coordinates
 *                    X = (p + xi)/2, Y = (q + eta)/2 (exact when the grid is uniform
 *                    in the coordinates in which the field is reconstructed)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_whitneyelement_new(WhitneyElement** self, int order, Grid_t* grid, int biquadratic);

/**
 * Destructor
 * @param self the element
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_whitneyelement_del(WhitneyElement** self);

/**
 * Get the order of the element
 * @param self the element
 * @param order 1 or 2 (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_whitneyelement_getOrder(WhitneyElement** self, int* order);

/**
 * Get the covariant components (a, b) of a field, a dxi + b deta in the
 * parametric coordinates of the given cells
 * @param self the element
 * @param grid the grid the data live on
 * @param numPoints number of points
 * @param cellIds cell of each point
 * @param xis first parametric coordinate of each point in its cell
 * @param etas second parametric coordinate of each point in its cell
 * @param data edge integrals, cell by cell (numCells*4, each cell's edges in the
 *             positive xi/eta direction) or on unique edges
 * @param placement MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA
 * @param a dxi components (output)
 * @param b deta components (output)
 * @return number of points whose cell is invalid or not covered by the element
 *         (their components are set to zero)
 */
LIBRARY_API
int mnt_whitneyelement_getCovariantComponents(WhitneyElement** self, Grid_t* grid, std::size_t numPoints,
                                              const long long cellIds[], const double xis[],
                                              const double etas[], const double data[], int placement,
                                              double a[], double b[]);

#endif // MNT_WHITNEY_ELEMENT
