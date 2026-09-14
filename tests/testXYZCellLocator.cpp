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

int main() {
    testXConstQuad();
    testTiltedQuad();
    std::cout << "SUCCESS\n";
    return 0;
}
