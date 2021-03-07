#include <mntPolysegmentIter.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vmtCellLocator.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>


void testFold() {
    // test folding at the pole

    // 2 cell mesh
    double v0[] = {  0., -90., 0.};
    double v1[] = {180., -90., 0.};
    double v2[] = {180.,  90., 0.};
    double v3[] = {  0.,  90., 0.};
    double v4[] = {180., -90., 0.};
    double v5[] = {360., -90., 0.};
    double v6[] = {360.,  90., 0.};
    double v7[] = {180.,  90., 0.};
    vtkPoints* points = vtkPoints::New();
    points->SetDataTypeToDouble();
    points->InsertNextPoint(v0);
    points->InsertNextPoint(v1);
    points->InsertNextPoint(v2);
    points->InsertNextPoint(v3);
    points->InsertNextPoint(v4);
    points->InsertNextPoint(v5);
    points->InsertNextPoint(v6);
    points->InsertNextPoint(v7);

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    size_t ncells = 2;
    grid->Allocate(ncells, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (size_t j = 0; j < ncells; ++j) {
       for (vtkIdType i = 0; i < 4; ++i) {
          ptIds->SetId(i, j*ncells + i);
       }
       grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    
    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {160.,  80.0, 0.};
    const double p1[] = {160., 100.0, 0.};
    const double p2[] = {200., 100.0, 0.};
    double xPeriod = 360.0;

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2)};
    assert(pts.size() >= 2);
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];

        PolysegmentIter psi(grid, loc, &pa[0], &pb[0], xPeriod);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "testFold: seg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1]
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1]
                                   << '\n';
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    ptIds->Delete();
    loc->Delete();
    points->Delete();
    grid->Delete();
}

void testFold2() {
    // test folding at the pole

    // 2 cell mesh
    double v0[] = {  0., -90., 0.};
    double v1[] = {180., -90., 0.};
    double v2[] = {180.,  90., 0.};
    double v3[] = {  0.,  90., 0.};
    double v4[] = {180., -90., 0.};
    double v5[] = {360., -90., 0.};
    double v6[] = {360.,  90., 0.};
    double v7[] = {180.,  90., 0.};
    vtkPoints* points = vtkPoints::New();
    points->SetDataTypeToDouble();
    points->InsertNextPoint(v0);
    points->InsertNextPoint(v1);
    points->InsertNextPoint(v2);
    points->InsertNextPoint(v3);
    points->InsertNextPoint(v4);
    points->InsertNextPoint(v5);
    points->InsertNextPoint(v6);
    points->InsertNextPoint(v7);

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    size_t ncells = 2;
    grid->Allocate(ncells, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (size_t j = 0; j < ncells; ++j) {
       for (vtkIdType i = 0; i < 4; ++i) {
          ptIds->SetId(i, j*ncells + i);
       }
       grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    
    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {-10.,  -90., 0.};
    const double p1[] = {-10.,   90., 0.};
    const double p2[] = {-10.,  100., 0.};
    const double p3[] = {370.,  100., 0.};
    const double p4[] = {370., -100., 0.};
    const double p5[] = {-10., -100.0, 0.};
    const double p6[] = {-10.,  -90.0, 0.};
    double xPeriod = 360.0;

    std::vector<Vec3> pts{Vec3(p0), Vec3(p1), Vec3(p2), Vec3(p3), Vec3(p4), Vec3(p5), Vec3(p6)};
    assert(pts.size() >= 2);
    for (size_t iInterval = 0; iInterval < pts.size() - 1; ++iInterval) {

        Vec3& pa = pts[iInterval + 0];
        Vec3& pb = pts[iInterval + 1];
        PolysegmentIter psi(grid, loc, &pa[0], &pb[0], xPeriod);

        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "testFold2: seg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1]
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1]
                                   << '\n';
            psi.next();
        }
        double tTotal = psi.getIntegratedParamCoord();
        double error = tTotal - 1.0;
        std::cout << "testFold: total t = " << tTotal <<  " error = " << error << '\n';
        assert(std::abs(error) < 1.e-10);
    }

    ptIds->Delete();
    loc->Delete();
    points->Delete();
    grid->Delete();
}


int main(int argc, char** argv) {

    //testFold2();
    testFold();

    return 0;
}
