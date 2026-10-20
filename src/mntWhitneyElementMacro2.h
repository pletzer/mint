#include <vector>
#include <mntWhitneyElement.h>

#ifndef MNT_WHITNEY_ELEMENT_MACRO2
#define MNT_WHITNEY_ELEMENT_MACRO2

/**
 * Second order edge element built on 2x2 blocks of grid cells (Rapetti-Bossavit
 * "small edge" degrees of freedom).
 *
 * The grid is tiled into blocks of 2x2 cells. In macro coordinates
 * (X, Y) in [0, 1]^2 of a block, the field is a Nedelec (first kind) degree 2
 * edge element, a dX + b dY with a in Q_{1,2} and b in Q_{2,1}. Its 12 degrees
 * of freedom are the integrals over the 12 edges of the block's grid cells, so
 * the element uses the same edge data as WhitneyElementOrder1:
 *
 *    + ---4---- + ---5---- +         Y
 *    |          |          |         ^
 *    7  (0,1)   9  (1,1)  11         |
 *    |          |          |         +--> X
 *    + ---2---- + ---3---- +
 *    |          |          |
 *    6  (0,0)   8  (1,0)  10
 *    |          |          |
 *    + ---0---- + ---1---- +
 *
 * (dof numbering: X-directed edges 0-5, row by row; Y-directed edges 6-11,
 * column by column, see getDofEdge.) The basis is
 *
 *    w_{2j + i}     = phi_i(X) L_j(Y) dX,   i = 0, 1   j = 0, 1, 2
 *    w_{6 + 2i + j} = L_i(X) phi_j(Y) dY,
 *
 * with L_j the quadratic Lagrange polynomials at 0, 1/2, 1 and
 * phi_0 = 3 - 4t, phi_1 = 4t - 1 the duals of the half-interval integrals.
 *
 * A grid cell at position (p, q) in its block has X = (p + xi)/2,
 * Y = (q + eta)/2, so the covariant components returned in the cell's own
 * parametric coordinates are half the macro components. The geometry is the
 * grid's own (piecewise bilinear) cell mapping.
 *
 * The four cells of a block must share their parametric frame (cell (1,0)
 * starts where cell (0,0)'s xi edge ends, etc.) and neighbouring blocks must
 * share their macro edges. build() fails on grids that cannot be tiled this
 * way, e.g. an odd number of cells along a direction.
 */
class WhitneyElementMacro2 : public WhitneyElement {

public:

    WhitneyElementMacro2() {}

    /**
     * Tile the grid into 2x2 blocks of cells
     * @param grid grid, must have UGRID face-node connectivity
     * @return error code (0 = OK), > 0 is the number of cells that could not be tiled
     */
    int build(Grid_t* grid);

    int getOrder() const override { return 2; }

    int getNumDofs() const override { return 12; }

    void evalBasis(vtkIdType cellId, double xi, double eta, double a[], double b[]) const override;

    void integrateBasis(vtkIdType cellId, const Vec3& xia, const Vec3& xib, double w[]) const override;

    int getCounterclockwiseSign(int dofIndex) const override;

    int getDofDataIndex(Grid_t* grid, vtkIdType cellId, int dofIndex, int placement,
                        std::size_t* index, int* sign) const override;

    /**
     * @return number of 2x2 blocks
     */
    std::size_t getNumBlocks() const { return this->blockCells.size() / 4; }

    /**
     * Get the cells of a block
     * @param blockId block Id
     * @param cellIds cells at positions (0,0), (1,0), (1,1), (0,1) (output)
     */
    void getBlockCells(std::size_t blockId, vtkIdType cellIds[]) const;

    /**
     * Get the block of a cell
     * @param cellId cell Id
     * @param blockId block Id (output)
     * @param p position of the cell in the block along X, 0 or 1 (output)
     * @param q position of the cell in the block along Y, 0 or 1 (output)
     * @return error code (0 = OK, 1 if the cell is not in any block)
     */
    int getCellBlock(vtkIdType cellId, std::size_t* blockId, int* p, int* q) const;

    /**
     * Get the grid cell edge of a degree of freedom
     * @param dofIndex local degree of freedom index
     * @param p position along X of the cell holding the edge (output)
     * @param q position along Y of the cell holding the edge (output)
     * @param edgeIndex edge index (0-3) in that cell (output)
     */
    static void getDofEdge(int dofIndex, int* p, int* q, int* edgeIndex);

    /**
     * Evaluate the basis functions in macro coordinates
     * @param x macro coordinate X in [0, 1]
     * @param y macro coordinate Y in [0, 1]
     * @param a dX components, array of size 12 (output)
     * @param b dY components, array of size 12 (output)
     */
    static void evalMacroBasis(double x, double y, double a[], double b[]);

private:

    // flat array of size numBlocks*4, cells at positions (0,0), (1,0), (1,1), (0,1)
    std::vector<vtkIdType> blockCells;

    // block Id of each cell, -1 if the cell is not in a block
    std::vector<long long> cellBlock;

    // position of each cell in its block, index into the 4 positions above
    std::vector<int> cellPos;
};

#endif // MNT_WHITNEY_ELEMENT_MACRO2
