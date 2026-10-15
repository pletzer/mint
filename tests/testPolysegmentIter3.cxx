#include <limits> // required by vtkUnstructuredGrid
#include <mntPolysegmentIter.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vmtLonLatCellLocator.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>

/**
 * Create uniform VTK grid objects
 * @param nx number of x cells
 * @param ny number of y cells
 * @param grid vtkUnstructured grid object (will be built)
 * @param points vtkPoints object (will be built)
 * @param coords vtkDoubleArray object (will be built)
 */
void createUniformGrid(int nx, int ny, 
                       vtkUnstructuredGrid* grid, vtkPoints* points, vtkDoubleArray* coords) {

    double point[3];

    int numCells = nx * ny;
    int numPoints = 4 * numCells;
    double dx = 360. / (double) nx;
    double dy = 180. / (double) ny;

    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(numPoints);

    int k = 0;
    for (int i = 0; i < nx; ++i) {
        double x0 = i*dx;
        double x1 = x0 + dx;
        for (int j = 0; j < ny; ++j) {
            double y0 = -90.0 + j*dy;
            double y1 = y0 + dy;
            point[0] = x0; point[1] = y0; point[2] = 0.0;
            coords->SetTuple(k*4 + 0, point);
            point[0] = x1; point[1] = y0; point[2] = 0.0;
            coords->SetTuple(k*4 + 1, point);
            point[0] = x1; point[1] = y1; point[2] = 0.0;
            coords->SetTuple(k*4 + 2, point);
            point[0] = x0; point[1] = y1; point[2] = 0.0;
            coords->SetTuple(k*4 + 3, point);
            k++;
        }
    }

    points->SetData(coords);

    grid->SetPoints(points);
    grid->Allocate(numCells, 1);
    const vtkIdType nptsPerCell = 4;
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(nptsPerCell);
    for (size_t iCell = 0; iCell < numCells; ++iCell) {
        for (size_t i = 0; i < nptsPerCell; ++i) {
            ptIds->SetId(i, nptsPerCell*iCell + i);
        }
        grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    ptIds->Delete();

}

void testSimple() {
    
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(36, 18, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    const double xPeriod = 360.0;
    loc->setPeriodicityLengthX(xPeriod);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {0.0,  0.0, 0.};
    const double p1[] = {10.0, 0.0, 0.};

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1)};
    assert(pts.size() >= 2);

    printf("testSimple\n");
    printf("intrvl    sgmt       cell       ta        tb       xia              xib\n");
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];

        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            // the parametric coords must be inside the cell
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            printf("%4zu     %4zu     %6lld   %.5f   %.5f  %.5f,%.5f  %.5f,%.5f\n", \
                    iInterval, i, cellId, ta, tb, xia[0], xia[1], xib[0], xib[1]);
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

void testFold() {
    // test folding at the pole

    
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(10, 5, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    const double xPeriod = 360.0;
    loc->setPeriodicityLengthX(xPeriod);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {160.,  80.0, 0.};
    const double p1[] = {160., 100.0, 0.};
    const double p2[] = {200., 100.0, 0.};

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2)};
    assert(pts.size() >= 2);

    printf("testFold\n");
    printf("intrvl    sgmt       cell       ta        tb       xia              xib\n");
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];

        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            // the parametric coords must be inside the cell
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            printf("%4zu     %4zu     %6lld   %.5f   %.5f  %.5f,%.5f  %.5f,%.5f\n", \
                    iInterval, i, cellId, ta, tb, xia[0], xia[1], xib[0], xib[1]);
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

void testFold2() {
    // test folding at the pole
    
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(10, 5, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    double xPeriod = 360.0;
    loc->setPeriodicityLengthX(xPeriod);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = { 10.,   80., 0.};
    const double p1[] = { 10.,  100., 0.};
    const double p2[] = {350.,  100., 0.};
    const double p3[] = {350.,   80., 0.};
    const double p4[] = {350.,  -80., 0.};
    const double p5[] = {350., -100., 0.};
    const double p6[] = { 10., -100., 0.};
    const double p7[] = { 10.,  -80., 0.};

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2), Vec3(p3), Vec3(p4), Vec3(p5), Vec3(p6), Vec3(p7)};
    assert(pts.size() >= 2);

    printf("testFold2\n");
    printf("intrvl    sgmt       cell       ta        tb       xia              xib\n");
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];
        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            // the parametric coords must be inside the cell
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            printf("%4zu     %4zu     %6lld   %.5f   %.5f  %.5f,%.5f  %.5f,%.5f\n", \
                    iInterval, i, cellId, ta, tb, xia[0], xia[1], xib[0], xib[1]);
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold2: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}


void testFold3() {
    // just outside of the domain
    
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(10, 5, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    double xPeriod = 360.0;
    loc->setPeriodicityLengthX(xPeriod);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = { -0.1, -90.1, 0.};
    const double p1[] = { 360.1, -90.1, 0.};
    const double p2[] = { 360.1, +90.1, 0.};
    const double p3[] = { -0.1,   +90.1, 0.};

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2), Vec3(p3)};
    assert(pts.size() >= 2);

    printf("testFold3\n");
    printf("intrvl    sgmt       cell       ta        tb       xia              xib\n");
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];
        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            // the parametric coords must be inside the cell
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            printf("%4zu     %4zu     %6lld   %.5f   %.5f  %.5f,%.5f  %.5f,%.5f\n", \
                    iInterval, i, cellId, ta, tb, xia[0], xia[1], xib[0], xib[1]);
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold3: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}


void testFold4() {
    // along the domain
    
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(10, 5, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    double xPeriod = 360.0;
    loc->setPeriodicityLengthX(xPeriod);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = { -0., -90., 0.};
    const double p1[] = { 360., -90., 0.};
    const double p2[] = { 360., +90., 0.};
    const double p3[] = { -0.,   +90., 0.};

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2), Vec3(p3)};
    assert(pts.size() >= 2);

    printf("testFold4\n");
    printf("intrvl    sgmt       cell       ta        tb       xia              xib\n");
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];
        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            // the parametric coords must be inside the cell
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            printf("%4zu     %4zu     %6lld   %.5f   %.5f  %.5f,%.5f  %.5f,%.5f\n", \
                    iInterval, i, cellId, ta, tb, xia[0], xia[1], xib[0], xib[1]);
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold4: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}


void testFold5() {
    // segments starting on the other side of the equator from the pole they cross. Both
    // end points must be folded about the same (crossed) pole, otherwise the folded line
    // is not the image of the original segment and pieces land in the wrong cells

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();

    createUniformGrid(36, 18, grid, points, coords);

    vmtLonLatCellLocator* loc = vmtLonLatCellLocator::New();
    loc->setPeriodicityLengthX(360.0);
    loc->enableFolding();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double segs[][4] = {
        // lon0, lat0, lon1, lat1
        {30., -80.,   30.,  100.},
        {30.,  -1.,   30.,  100.},
        {30., -1.e-6, 30.,  100.},
        {30.,   0.,   30.,  100.},
        {30.,  80.,   30., -100.},
        {30.,   1.,   30., -100.},
        {30., -80.,   50.,  100.},
        {200., 20.,  190., -110.},
        {30., -95.,   30.,   95.},  // beyond both poles
        {30.,  95.,   50.,  -95.},
    };

    printf("testFold5\n");
    for (const auto& s : segs) {

        Vec3 pa, pb;
        pa[0] = s[0]; pa[1] = s[1]; pa[2] = 0.;
        pb[0] = s[2]; pb[1] = s[3]; pb[2] = 0.;
        PolysegmentIter psi(grid, loc, &pa[0], &pb[0]);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            for (size_t k = 0; k < 2; ++k) {
                assert(xia[k] > -1.e-10 && xia[k] < 1. + 1.e-10);
                assert(xib[k] > -1.e-10 && xib[k] < 1. + 1.e-10);
            }

            // the mid point of this piece of the original segment, brought back into
            // the domain, must lie in the cell
            double t = 0.5*(psi.getBegLineParamCoord() + psi.getEndLineParamCoord());
            double lon = pa[0] + t*(pb[0] - pa[0]);
            double lat = pa[1] + t*(pb[1] - pa[1]);
            if (lat > 90.) {
                lat = 180. - lat;
                lon += 180.;
            }
            else if (lat < -90.) {
                lat = -180. - lat;
                lon += 180.;
            }
            lon -= 360.*std::floor(lon/360.);
            double* bounds = grid->GetCell(cellId)->GetBounds();
            const double tol = 1.e-8;
            if (lon < bounds[0] - tol || lon > bounds[1] + tol ||
                lat < bounds[2] - tol || lat > bounds[3] + tol) {
                printf("testFold5: segment (%g,%g)->(%g,%g) t = %g: mid point (%g,%g) not in cell %lld [%g,%g]x[%g,%g]\n",
                       s[0], s[1], s[2], s[3], t, lon, lat, cellId,
                       bounds[0], bounds[1], bounds[2], bounds[3]);
                assert(false);
            }
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        printf("testFold5: (%g,%g)->(%g,%g) total t = %.12f error = %g\n", s[0], s[1], s[2], s[3], tTotal, error);
        assert(std::abs(error) < 1.e-10);
    }

    loc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();
}

int main(int argc, char** argv) {

    testSimple();
    testFold5();
    testFold4();
    testFold3();
    testFold2();
    testFold();

    return 0;
}
