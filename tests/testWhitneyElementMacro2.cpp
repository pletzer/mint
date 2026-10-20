#include <mntWhitneyElementMacro2.h>
#include <mntGlobal.h>
#include <mntGrid.h>
#include <mntVectorInterp.h>
#include <mntPolylineIntegral.h>
#include <mntRegridEdges.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <cstring>

// affine map from grid index space (s, t) to physical (x, y), so that a
// potential quadratic in (x, y) is quadratic in the macro coordinates
void affineMap(double s, double t, double* x, double* y) {
    *x = 0.2 + 1.0*s + 0.25*t;
    *y = -0.1 + 0.15*s + 0.8*t;
}

double potential(double x, double y) {
    return 0.3*x*x - 0.7*x*y + 0.2*y*y + 0.5*x - 1.1*y;
}

void gradPotential(double x, double y, double g[]) {
    g[0] = 0.6*x - 0.7*y + 0.5;
    g[1] = -0.7*x + 0.4*y - 1.1;
}

// node pairs (local node indices) of the 4 cell edges, positive xi/eta direction
const int EDGE_NODES[] = {0, 1,   1, 2,   3, 2,   0, 3};

/**
 * Build a structured nx x ny UGRID grid of the box [s0, s1] x [t0, t1] in index
 * space, mapped by affineMap
 */
Grid_t* createGrid(std::size_t nx, std::size_t ny,
                   double s0 = 0., double s1 = -1., double t0 = 0., double t1 = -1.) {
    if (s1 < s0) s1 = (double) nx;
    if (t1 < t0) t1 = (double) ny;
    std::vector<double> xyz;
    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            double x, y;
            affineMap(s0 + (s1 - s0)*i/nx, t0 + (t1 - t0)*j/ny, &x, &y);
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(0.0);
        }
    }
    auto node = [nx](std::size_t i, std::size_t j) { return j*(nx + 1) + i; };
    std::vector<std::size_t> f2n, e2n;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            f2n.push_back(node(i, j));
            f2n.push_back(node(i + 1, j));
            f2n.push_back(node(i + 1, j + 1));
            f2n.push_back(node(i, j + 1));
        }
    }
    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            e2n.push_back(node(i, j));
            e2n.push_back(node(i + 1, j));
        }
    }
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            e2n.push_back(node(i, j));
            e2n.push_back(node(i, j + 1));
        }
    }
    Grid_t* grid = NULL;
    mnt_grid_new(&grid);
    mnt_grid_setFlags(&grid, 0, 0, 0);
    int ier = mnt_grid_loadFromUgrid2DData(&grid, nx*ny, e2n.size()/2, xyz.size()/3,
                                           &xyz[0], &f2n[0], &e2n[0]);
    assert(ier == 0);
    return grid;
}

/**
 * Edge integrals of the gradient of the potential, cell by cell and on unique edges
 */
void getEdgeData(Grid_t* grid, std::vector<double>& cellByCell, std::vector<double>& uniqueEdge) {
    std::size_t numCells, numEdges;
    mnt_grid_getNumberOfCells(&grid, &numCells);
    mnt_grid_getNumberOfEdges(&grid, &numEdges);
    const double* xyz = grid->verts;
    cellByCell.resize(numCells*4);
    for (std::size_t c = 0; c < numCells; ++c) {
        for (int e = 0; e < 4; ++e) {
            const double* p0 = &xyz[3*(4*c + EDGE_NODES[2*e])];
            const double* p1 = &xyz[3*(4*c + EDGE_NODES[2*e + 1])];
            cellByCell[4*c + e] = potential(p1[0], p1[1]) - potential(p0[0], p0[1]);
        }
    }
    uniqueEdge.assign(numEdges, 0.0);
    for (std::size_t c = 0; c < numCells; ++c) {
        for (int e = 0; e < 4; ++e) {
            std::size_t edgeId;
            int edgeSign;
            mnt_grid_getEdgeId(&grid, c, e, &edgeId, &edgeSign);
            uniqueEdge[edgeId] = edgeSign * cellByCell[4*c + e];
        }
    }
}

/**
 * The grid cell edge (p, q, edgeIndex) of a block, made unique: an edge shared
 * by two cells of the block is referred to by the cell at the lower position
 */
void canonicalEdge(int* p, int* q, int* e) {
    if (*e == 0 && *q == 1) { *q = 0; *e = 2; }
    if (*e == 3 && *p == 1) { *p = 0; *e = 1; }
}

void testTiling() {
    Grid_t* grid = createGrid(6, 4);
    WhitneyElementMacro2 elem;
    assert(elem.build(grid) == 0);
    assert(elem.getNumBlocks() == 6);
    for (vtkIdType c = 0; c < 24; ++c) {
        std::size_t blockId;
        int p, q;
        assert(elem.getCellBlock(c, &blockId, &p, &q) == 0);
        vtkIdType cells[4];
        elem.getBlockCells(blockId, cells);
        assert(cells[(q == 0)? p: 3 - p] == c);
    }
    mnt_grid_del(&grid);

    // an odd number of cells cannot be tiled
    grid = createGrid(5, 4);
    WhitneyElementMacro2 elem2;
    assert(elem2.build(grid) > 0);
    mnt_grid_del(&grid);

    std::cout << "testTiling: OK\n";
}

void testTilingMeshes() {
    // cubed-sphere meshes tile across the panel edges, the lon-lat mesh across the poles
    const char* files[] = {"${CMAKE_SOURCE_DIR}/data/cs_16.nc$physics",
                           "${CMAKE_SOURCE_DIR}/data/lfric_diag_wind.nc$Mesh2d",
                           "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc$latlon"};
    const std::size_t numBlocks[] = {384, 384, 1250};
    for (int i = 0; i < 3; ++i) {
        Grid_t* grid = NULL;
        mnt_grid_new(&grid);
        int ier = mnt_grid_loadFromUgrid2DFile(&grid, files[i]);
        assert(ier == 0);
        WhitneyElementMacro2 elem;
        assert(elem.build(grid) == 0);
        assert(elem.getNumBlocks() == numBlocks[i]);
        mnt_grid_del(&grid);
    }
    std::cout << "testTilingMeshes: OK\n";
}

void testDuality() {
    // the integral of basis function k over the cell edge of degree of freedom j is delta_kj,
    // seen from either cell sharing that edge
    const double ends[] = {0., 0., 0.,  1., 0., 0.,  1., 1., 0.,  0., 1., 0.};
    Grid_t* grid = createGrid(2, 2);
    WhitneyElementMacro2 elem;
    assert(elem.build(grid) == 0);
    double w[12];
    for (vtkIdType c = 0; c < 4; ++c) {
        std::size_t blockId;
        int p, q;
        elem.getCellBlock(c, &blockId, &p, &q);
        for (int e = 0; e < 4; ++e) {
            Vec3 xia, xib;
            for (int d = 0; d < 3; ++d) {
                xia[d] = ends[3*EDGE_NODES[2*e] + d];
                xib[d] = ends[3*EDGE_NODES[2*e + 1] + d];
            }
            elem.integrateBasis(c, xia, xib, w);
            int pc = p, qc = q, ec = e;
            canonicalEdge(&pc, &qc, &ec);
            for (int k = 0; k < 12; ++k) {
                int pk, qk, ek;
                WhitneyElementMacro2::getDofEdge(k, &pk, &qk, &ek);
                double expected = (pk == pc && qk == qc && ek == ec)? 1.0: 0.0;
                assert(std::abs(w[k] - expected) < 1.e-14);
            }
        }
    }
    mnt_grid_del(&grid);
    std::cout << "testDuality: OK\n";
}

void testIntegrateBasis() {
    // exact segment integrals agree with quadrature of the basis functions
    Grid_t* grid = createGrid(2, 2);
    WhitneyElementMacro2 elem;
    assert(elem.build(grid) == 0);
    std::mt19937 gen(123);
    std::uniform_real_distribution<double> u(0., 1.);
    double w[12], wq[12], a[12], b[12];
    for (vtkIdType c = 0; c < 4; ++c) {
        for (int iseg = 0; iseg < 5; ++iseg) {
            Vec3 xia, xib;
            xia[0] = u(gen); xia[1] = u(gen); xia[2] = 0.;
            xib[0] = u(gen); xib[1] = u(gen); xib[2] = 0.;
            elem.integrateBasis(c, xia, xib, w);
            for (int k = 0; k < 12; ++k) wq[k] = 0.;
            const int n = 2000;
            for (int m = 0; m < n; ++m) {
                double t = (m + 0.5)/n;
                elem.evalBasis(c, xia[0] + t*(xib[0] - xia[0]), xia[1] + t*(xib[1] - xia[1]), a, b);
                for (int k = 0; k < 12; ++k) {
                    wq[k] += (a[k]*(xib[0] - xia[0]) + b[k]*(xib[1] - xia[1]))/n;
                }
            }
            for (int k = 0; k < 12; ++k) {
                assert(std::abs(w[k] - wq[k]) < 1.e-6);
            }
        }
    }
    mnt_grid_del(&grid);
    std::cout << "testIntegrateBasis: OK\n";
}

void testTangentialContinuity() {
    // inside a block, the tangential component is continuous across the cell edges
    Grid_t* grid = createGrid(2, 2);
    WhitneyElementMacro2 elem;
    assert(elem.build(grid) == 0);
    const double dofs[] = {1.3, -0.7, 2.1, 0.4, 0.9, -1.1, 0.2, 0.5, -0.3, 1.7, 0.8, -2.0};
    vtkIdType cells[4];
    elem.getBlockCells(0, cells);
    for (double t = 0.; t <= 1.0; t += 0.125) {
        double a0, b0, a1, b1;
        // (0,0) | (1,0): deta component
        elem.getCovariantComponents(cells[0], 1.0, t, dofs, &a0, &b0);
        elem.getCovariantComponents(cells[1], 0.0, t, dofs, &a1, &b1);
        assert(std::abs(b0 - b1) < 1.e-14);
        // (0,0) / (0,1): dxi component
        elem.getCovariantComponents(cells[0], t, 1.0, dofs, &a0, &b0);
        elem.getCovariantComponents(cells[3], t, 0.0, dofs, &a1, &b1);
        assert(std::abs(a0 - a1) < 1.e-14);
    }
    mnt_grid_del(&grid);
    std::cout << "testTangentialContinuity: OK\n";
}

/**
 * Gradient of the quadratic potential at random points through VectorInterp,
 * @return max error of the edge vectors
 */
double vectorInterpError(int order, int placement, bool faceVectors) {

    Grid_t* grid = createGrid(6, 4);
    std::vector<double> cellByCell, uniqueEdge;
    getEdgeData(grid, cellByCell, uniqueEdge);

    VectorInterp_t* vp = NULL;
    mnt_vectorinterp_new(&vp);
    mnt_vectorinterp_setGrid(&vp, grid);
    int ier = mnt_vectorinterp_setOrder(&vp, order);
    assert(ier == 0);
    int order2;
    mnt_vectorinterp_getOrder(&vp, &order2);
    assert(order2 == order);
    mnt_vectorinterp_buildLocator(&vp, 8, 0., 0, 0);

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> us(0.01, 5.99), ut(0.01, 3.99);
    const std::size_t n = 200;
    std::vector<double> targets(3*n), vectors(3*n);
    for (std::size_t i = 0; i < n; ++i) {
        affineMap(us(gen), ut(gen), &targets[3*i], &targets[3*i + 1]);
        targets[3*i + 2] = 0.;
    }
    ier = mnt_vectorinterp_findPoints(&vp, n, &targets[0], 1.e-10);
    assert(ier == 0);
    const double* data = (placement == MNT_CELL_BY_CELL_DATA)? &cellByCell[0]: &uniqueEdge[0];
    if (faceVectors) {
        ier = mnt_vectorinterp_getFaceVectors(&vp, data, placement, &vectors[0]);
    }
    else {
        ier = mnt_vectorinterp_getEdgeVectors(&vp, data, placement, &vectors[0]);
    }
    assert(ier == 0);

    double error = 0.;
    for (std::size_t i = 0; i < n; ++i) {
        double g[2];
        gradPotential(targets[3*i], targets[3*i + 1], g);
        // the face vector is the edge vector rotated by -90 degrees
        double ex = faceVectors? g[1]: g[0];
        double ey = faceVectors? -g[0]: g[1];
        error = std::max(error, std::abs(vectors[3*i] - ex) + std::abs(vectors[3*i + 1] - ey));
    }
    mnt_vectorinterp_del(&vp);
    mnt_grid_del(&grid);
    return error;
}

void testVectorInterp() {
    for (int placement : {MNT_CELL_BY_CELL_DATA, MNT_UNIQUE_EDGE_DATA}) {
        for (bool faceVectors : {false, true}) {
            double error1 = vectorInterpError(1, placement, faceVectors);
            double error2 = vectorInterpError(2, placement, faceVectors);
            std::cout << "testVectorInterp: placement " << placement << " face " << faceVectors
                      << " order 1 error " << error1 << " macro error " << error2 << '\n';
            // the macro element reproduces the gradient of a quadratic exactly, order 1 does not
            assert(error2 < 1.e-10);
            assert(error1 > 1.e-3);
        }
    }
    std::cout << "testVectorInterp: OK\n";
}

/**
 * Line integral of the gradient of the potential along a polyline
 * @return error
 */
double polylineIntegralError(int order, int placement) {

    Grid_t* grid = createGrid(6, 4);
    std::vector<double> cellByCell, uniqueEdge;
    getEdgeData(grid, cellByCell, uniqueEdge);

    PolylineIntegral_t* pli = NULL;
    mnt_polylineintegral_new(&pli);
    mnt_polylineintegral_setGrid(&pli, grid);
    int ier = mnt_polylineintegral_setOrder(&pli, order);
    assert(ier == 0);
    mnt_polylineintegral_buildLocator(&pli, 8, 0., 0, 0);

    const double st[] = {0.3, 0.2,   2.7, 3.1,   4.9, 1.3,   5.6, 3.7};
    std::vector<double> xyz(12, 0.0);
    for (int i = 0; i < 4; ++i) {
        affineMap(st[2*i], st[2*i + 1], &xyz[3*i], &xyz[3*i + 1]);
    }
    ier = mnt_polylineintegral_computeWeights(&pli, 4, &xyz[0], 0);
    assert(ier == 0);
    const double* data = (placement == MNT_CELL_BY_CELL_DATA)? &cellByCell[0]: &uniqueEdge[0];
    double result;
    ier = mnt_polylineintegral_getIntegral(&pli, data, placement, &result);
    assert(ier == 0);
    double exact = potential(xyz[9], xyz[10]) - potential(xyz[0], xyz[1]);

    mnt_polylineintegral_del(&pli);
    mnt_grid_del(&grid);
    return std::abs(result - exact);
}

void testPolylineIntegral() {
    for (int placement : {MNT_CELL_BY_CELL_DATA, MNT_UNIQUE_EDGE_DATA}) {
        double error1 = polylineIntegralError(1, placement);
        double error2 = polylineIntegralError(2, placement);
        std::cout << "testPolylineIntegral: placement " << placement
                  << " order 1 error " << error1 << " macro error " << error2 << '\n';
        assert(error2 < 1.e-10);
        assert(error1 > 1.e-4);
    }
    std::cout << "testPolylineIntegral: OK\n";
}

/**
 * Regrid the edge integrals of the gradient of the potential onto a skewed grid
 * @return max error of the destination edge integrals
 */
double regridError(int order, const char* weightsFile = NULL) {

    Grid_t* srcGrid = createGrid(6, 4);
    Grid_t* dstGrid = createGrid(3, 3, 0.4, 5.3, 0.3, 3.6);
    std::vector<double> srcCellByCell, srcUniqueEdge, dstCellByCell, dstUniqueEdge;
    getEdgeData(srcGrid, srcCellByCell, srcUniqueEdge);
    getEdgeData(dstGrid, dstCellByCell, dstUniqueEdge);

    RegridEdges_t* rg = NULL;
    mnt_regridedges_new(&rg);
    mnt_regridedges_setSrcGrid(&rg, srcGrid);
    mnt_regridedges_setDstGrid(&rg, dstGrid);
    int ier = mnt_regridedges_setOrder(&rg, order);
    assert(ier == 0);
    if (weightsFile) {
        // weights computed earlier with the same order
        ier = mnt_regridedges_loadWeights(&rg, weightsFile, (int) strlen(weightsFile));
        assert(ier == 0);
    }
    else {
        mnt_regridedges_buildLocator(&rg, 8, 0., 0, 0);
        ier = mnt_regridedges_computeWeights(&rg, 0);
        assert(ier == 0);
    }
    std::vector<double> result(dstUniqueEdge.size());
    ier = mnt_regridedges_apply(&rg, &srcUniqueEdge[0], &result[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    double error = 0.;
    for (std::size_t i = 0; i < result.size(); ++i) {
        error = std::max(error, std::abs(result[i] - dstUniqueEdge[i]));
    }
    mnt_regridedges_del(&rg);
    mnt_grid_del(&srcGrid);
    mnt_grid_del(&dstGrid);
    return error;
}

void testRegridEdges() {
    double error1 = regridError(1);
    double error2 = regridError(2);
    std::cout << "testRegridEdges: order 1 error " << error1 << " macro error " << error2 << '\n';
    assert(error2 < 1.e-10);
    assert(error1 > 1.e-4);
    std::cout << "testRegridEdges: OK\n";
}

void testSetOrderErrors() {

    int order;

    // invalid order
    VectorInterp_t* vp = NULL;
    mnt_vectorinterp_new(&vp);
    assert(mnt_vectorinterp_setOrder(&vp, 3) != 0);
    // order 2 needs a grid
    assert(mnt_vectorinterp_setOrder(&vp, 2) != 0);
    mnt_vectorinterp_getOrder(&vp, &order);
    assert(order == 1);

    // a grid that cannot be tiled into 2x2 blocks
    Grid_t* oddGrid = createGrid(5, 4);
    mnt_vectorinterp_setGrid(&vp, oddGrid);
    assert(mnt_vectorinterp_setOrder(&vp, 2) != 0);
    mnt_vectorinterp_getOrder(&vp, &order);
    assert(order == 1);

    // setting a new grid rebuilds the order 2 element, or falls back to order 1
    Grid_t* grid = createGrid(6, 4);
    assert(mnt_vectorinterp_setGrid(&vp, grid) == 0);
    assert(mnt_vectorinterp_setOrder(&vp, 2) == 0);
    Grid_t* grid2 = createGrid(4, 4);
    assert(mnt_vectorinterp_setGrid(&vp, grid2) == 0);
    mnt_vectorinterp_getOrder(&vp, &order);
    assert(order == 2);
    assert(mnt_vectorinterp_setGrid(&vp, oddGrid) != 0);
    mnt_vectorinterp_getOrder(&vp, &order);
    assert(order == 1);
    mnt_vectorinterp_del(&vp);

    // changing the order discards the polyline weights
    PolylineIntegral_t* pli = NULL;
    mnt_polylineintegral_new(&pli);
    mnt_polylineintegral_setGrid(&pli, grid);
    mnt_polylineintegral_buildLocator(&pli, 8, 0., 0, 0);
    double xyz[] = {0.5, 0.5, 0.,  4.0, 2.0, 0.};
    assert(mnt_polylineintegral_computeWeights(&pli, 2, xyz, 0) == 0);
    assert(pli->weights.size() > 0);
    assert(mnt_polylineintegral_setOrder(&pli, 2) == 0);
    assert(pli->weights.size() == 0);
    mnt_polylineintegral_getOrder(&pli, &order);
    assert(order == 2);
    mnt_polylineintegral_del(&pli);

    // order 2 weights files can only be loaded with order 2
    RegridEdges_t* rg = NULL;
    mnt_regridedges_new(&rg);
    Grid_t* dstGrid = createGrid(3, 3, 0.4, 5.3, 0.3, 3.6);
    mnt_regridedges_setSrcGrid(&rg, grid);
    mnt_regridedges_setDstGrid(&rg, dstGrid);
    assert(mnt_regridedges_setOrder(&rg, 2) == 0);
    mnt_regridedges_buildLocator(&rg, 8, 0., 0, 0);
    assert(mnt_regridedges_computeWeights(&rg, 0) == 0);
    const char* weightsFile = "testWhitneyElementMacro2_weights.nc";
    assert(mnt_regridedges_dumpWeights(&rg, weightsFile, (int) strlen(weightsFile)) == 0);
    assert(mnt_regridedges_setOrder(&rg, 1) == 0);
    assert(rg->weights.size() == 0);
    assert(mnt_regridedges_loadWeights(&rg, weightsFile, (int) strlen(weightsFile)) != 0);
    mnt_regridedges_del(&rg);
    // loading them with order 2 gives the same (exact) result
    assert(regridError(2, weightsFile) < 1.e-10);

    mnt_grid_del(&oddGrid);
    mnt_grid_del(&grid);
    mnt_grid_del(&grid2);
    mnt_grid_del(&dstGrid);
    std::cout << "testSetOrderErrors: OK\n";
}

int main(int argc, char** argv) {

    testTiling();
    testTilingMeshes();
    testDuality();
    testIntegrateBasis();
    testTangentialContinuity();
    testVectorInterp();
    testPolylineIntegral();
    testRegridEdges();
    testSetOrderErrors();

    return 0;
}
