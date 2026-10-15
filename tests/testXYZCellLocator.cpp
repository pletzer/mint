#include <limits> // required by vtkUnstructuredGrid
#include <iostream>
#include <cmath>
#include <vector>
#include <mntGrid.h>
#include <vmtXYZCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#undef NDEBUG // turn on asserts
#include <cassert>

/**
 * Direct C++ (not Python/ctypes-mediated) regression tests for
 * vmtXYZCellLocator, mirroring testCellLocator.cpp's coverage of
 * vmtLonLatCellLocator. Before this file, vmtXYZCellLocator's only test
 * coverage was via Python (mint/tests/test_cell_locator_xyz.py,
 * test_vector_interp_xyz.py) -- useful, but it only proves correctness
 * under whatever VTK version the Python environment happens to link
 * against, which is not necessarily the same VTK version ctest's own
 * build uses (see the VTK 9.1-vs-9.4 vtkStaticCellLocator::FindCell
 * overload mismatch this project already hit -- a pure compile error that
 * pytest could never have caught, but this file also gives ctest's build
 * a direct RUNTIME check of vmtXYZCellLocator under its own VTK version,
 * not just a compile check).
 */

/**
 * Build a one-cell VTK grid from 4 explicit (x, y, z) corners.
 */
void buildSingleCellGrid(const double p0[3], const double p1[3], const double p2[3], const double p3[3],
                         vtkUnstructuredGrid* grid, vtkPoints* points, vtkDoubleArray* coords) {
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(4);
    coords->SetTuple(0, p0);
    coords->SetTuple(1, p1);
    coords->SetTuple(2, p2);
    coords->SetTuple(3, p3);
    points->SetData(coords);

    grid->SetPoints(points);
    grid->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    grid->InsertNextCell(VTK_QUAD, ptIds);
    ptIds->Delete();
}

/**
 * A cell lying entirely in the x=0 plane -- the case that used to give the
 * OLD (x,y)-only crossDotZHat/isPointInQuad an exactly-zero Jacobian (see
 * test_vector_interp_xyz.py's test_x_const_quad_is_not_degenerate for the
 * Python-level version). Checks vmtXYZCellLocator finds an interior point
 * and correctly rejects one well outside the quad's own (y, z) footprint.
 */
void testXConstQuad() {
    double p0[3] = {0., 0., 0.};
    double p1[3] = {0., 1., 0.};
    double p2[3] = {0., 1., 1.};
    double p3[3] = {0., 0., 1.};

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();
    buildSingleCellGrid(p0, p1, p2, p3, grid, points, coords);

    vmtXYZCellLocator* loc = vmtXYZCellLocator::New();
    loc->SetDataSet(grid);
    loc->SetNumberOfCellsPerBucket(1);
    loc->BuildLocator();

    const double tol2 = 1.e-10;
    vtkGenericCell* cell = NULL;
    double pcoords[3];
    std::vector<double> weights(4);

    double target[3] = {0., 0.5, 0.5}; // cell centre
    vtkIdType cellId = loc->FindCell(target, tol2, cell, pcoords, &weights[0]);
    assert(cellId == 0);

    double outside[3] = {0., 5., 5.}; // same x=0 plane, well outside the unit square
    vtkIdType badId = loc->FindCell(outside, tol2, cell, pcoords, &weights[0]);
    assert(badId < 0);

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

/**
 * A cell whose plane's normal, (1,1,1)/sqrt(3), isn't aligned with any
 * single coordinate axis -- unlike testXConstQuad, which just swaps which
 * axis is "dropped", this checks the locator on a genuinely tilted cell
 * (see test_vector_interp_xyz.py's test_tilted_quad_arbitrary_normal).
 */
void testTiltedQuad() {
    double s3 = std::sqrt(3.0), s2 = std::sqrt(2.0);
    double n[3]  = {1./s3, 1./s3, 1./s3};
    double e1[3] = {1./s2, -1./s2, 0.};
    double e2[3] = {n[1]*e1[2] - n[2]*e1[1],
                    n[2]*e1[0] - n[0]*e1[2],
                    n[0]*e1[1] - n[1]*e1[0]};
    double base[3] = {1., 1., 1.};

    double p0[3], p1[3], p2[3], p3[3];
    for (int i = 0; i < 3; ++i) {
        p0[i] = base[i];
        p1[i] = base[i] + e1[i];
        p2[i] = base[i] + e1[i] + e2[i];
        p3[i] = base[i] + e2[i];
    }

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();
    buildSingleCellGrid(p0, p1, p2, p3, grid, points, coords);

    vmtXYZCellLocator* loc = vmtXYZCellLocator::New();
    loc->SetDataSet(grid);
    loc->SetNumberOfCellsPerBucket(1);
    loc->BuildLocator();

    // cell centre: bilinear position at (xi, eta) = (0.5, 0.5), i.e. the
    // plain average of the 4 corners
    double target[3];
    for (int i = 0; i < 3; ++i) {
        target[i] = 0.25*(p0[i] + p1[i] + p2[i] + p3[i]);
    }

    const double tol2 = 1.e-10;
    vtkGenericCell* cell = NULL;
    double pcoords[3];
    std::vector<double> weights(4);
    vtkIdType cellId = loc->FindCell(target, tol2, cell, pcoords, &weights[0]);
    assert(cellId == 0);
    assert(std::abs(pcoords[0] - 0.5) < 1.e-6);
    assert(std::abs(pcoords[1] - 0.5) < 1.e-6);

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

/**
 * Trace lines across a flat n x n grid of unit-square cells in the z=0 plane
 * with findIntersectionsWithLine. The segments must cover the whole line,
 * [0, 1] in t, without gaps: each segment must start where the previous one
 * ended, even though the next cell is looked up slightly past the boundary
 */
void testLineAcrossCells() {
    const int n = 10;
    const double h = 1.0 / n;

    vtkDoubleArray* coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(4*n*n);
    vtkIdType k = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            double x0 = i*h, x1 = (i + 1)*h, y0 = j*h, y1 = (j + 1)*h;
            double c[4][3] = {{x0, y0, 0.}, {x1, y0, 0.}, {x1, y1, 0.}, {x0, y1, 0.}};
            for (int v = 0; v < 4; ++v) {
                coords->SetTuple(k*4 + v, c[v]);
            }
            k++;
        }
    }
    vtkPoints* points = vtkPoints::New();
    points->SetData(coords);
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    grid->Allocate(n*n, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType iCell = 0; iCell < n*n; ++iCell) {
        for (vtkIdType v = 0; v < 4; ++v) {
            ptIds->SetId(v, iCell*4 + v);
        }
        grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    ptIds->Delete();

    vmtXYZCellLocator* loc = vmtXYZCellLocator::New();
    loc->SetDataSet(grid);
    loc->SetNumberOfCellsPerBucket(4);
    loc->BuildLocator();

    const double lines[][6] = {
        {0.05, 0.03, 0.,  0.95,    0.87, 0.},
        // ends just past a cell boundary, the last cell is entered within
        // less than the look-up nudge of the end of the line
        {0.05, 0.55, 0.,  0.90005, 0.55, 0.},
    };
    for (const auto& line : lines) {
        Vec3 pBeg(&line[0]), pEnd(&line[3]);
        std::vector< std::pair<vtkIdType, Vec4> > segs = loc->findIntersectionsWithLine(pBeg, pEnd);
        assert(segs.size() > 0);

        double totalT = 0.0;
        double tPrev = 0.0;
        for (const auto& s : segs) {
            double ta = s.second[0], tb = s.second[1];
            // no gap or overlap with the previous segment
            assert(std::abs(ta - tPrev) < 1.e-12);
            tPrev = tb;
            totalT += tb - ta;

            // the segment's mid point lies in its cell
            Vec3 pMid = pBeg + (pEnd - pBeg)*(0.5*(ta + tb));
            double* bounds = grid->GetCell(s.first)->GetBounds();
            assert(pMid[0] > bounds[0] - 1.e-6 && pMid[0] < bounds[1] + 1.e-6);
            assert(pMid[1] > bounds[2] - 1.e-6 && pMid[1] < bounds[3] + 1.e-6);
        }
        std::cout << "testLineAcrossCells: " << segs.size() << " segments, total t = "
                  << totalT << " error = " << totalT - 1.0 << '\n';
        assert(std::abs(totalT - 1.0) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

int main() {
    testXConstQuad();
    testTiltedQuad();
    testLineAcrossCells();
    std::cout << "SUCCESS\n";
    return 0;
}
