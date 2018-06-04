#include <mntPolysegmentIter.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellLocator.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>

void test1Cell() {
    // a one cell grid, built from scratch
    double v0[] = {0., 0., 0.};
    double v1[] = {1., 0., 0.};
    double v2[] = {1., 1., 0.};
    double v3[] = {0., 1., 0.};
    vtkPoints* points = vtkPoints::New();
    points->SetDataTypeToDouble();
    points->InsertNextPoint(v0);
    points->InsertNextPoint(v1);
    points->InsertNextPoint(v2);
    points->InsertNextPoint(v3);

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    grid->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    grid->InsertNextCell(VTK_QUAD, ptIds);
    

    vtkCellLocator* loc = vtkCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {0.2, 0.3, 0.};
    const double p1[] = {0.4, 0.5, 0.};
    PolysegmentIter psi(grid, loc, p0, p1);
    size_t numSegs = psi.getNumberOfSegments();
    psi.reset();
    for (size_t i = 0; i < numSegs; ++i) {
        vtkIdType cellId = psi.getCellId();
        const std::vector<double>& xia = psi.getBegCellParamCoord();
        const std::vector<double>& xib = psi.getEndCellParamCoord();
        double ta = psi.getBegLineParamCoord();
        double tb = psi.getEndLineParamCoord();
        double coeff = psi.getCoefficient();
        std::cout << "test1Cell: seg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                   << '\n';
        psi.next();
    }
    assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);

    ptIds->Delete();
    loc->Delete();
    grid->Delete();
    points->Delete();
}

void test2Cells() {
    // two cell grid, built from scratch
    double v0[] = {0., 0., 0.};
    double v1[] = {1., 0., 0.};
    double v2[] = {1., 1., 0.};
    double v3[] = {0., 1., 0.};
    double v4[] = {1., 0., 0.};
    double v5[] = {2., 0., 0.};
    double v6[] = {2., 1., 0.};
    double v7[] = {1., 1., 0.};
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
    
    vtkCellLocator* loc = vtkCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {0.2, 0.3, 0.};
    const double p1[] = {1.4, 0.5, 0.};
    PolysegmentIter psi(grid, loc, p0, p1);
    size_t numSegs = psi.getNumberOfSegments();
    psi.reset();
    for (size_t i = 0; i < numSegs; ++i) {
        vtkIdType cellId = psi.getCellId();
        const std::vector<double>& xia = psi.getBegCellParamCoord();
        const std::vector<double>& xib = psi.getEndCellParamCoord();
        double ta = psi.getBegLineParamCoord();
        double tb = psi.getEndLineParamCoord();
        double coeff = psi.getCoefficient();
        std::cout << "test2Cells: seg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                   << '\n';
        psi.next();
    }
    assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);

    ptIds->Delete();
    loc->Delete();
    grid->Delete();
    points->Delete();
}

void test2CellsEdge() {
    // two cell grid, built from scratch
    double v0[] = {0., 0., 0.};
    double v1[] = {1., 0., 0.};
    double v2[] = {1., 1., 0.};
    double v3[] = {0., 1., 0.};
    double v4[] = {1., 0., 0.};
    double v5[] = {2., 0., 0.};
    double v6[] = {2., 1., 0.};
    double v7[] = {1., 1., 0.};
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
    
    vtkCellLocator* loc = vtkCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {0.2, 0.0, 0.};
    const double p1[] = {1.4, 0.0, 0.};
    PolysegmentIter psi(grid, loc, p0, p1);
    size_t numSegs = psi.getNumberOfSegments();
    psi.reset();
    for (size_t i = 0; i < numSegs; ++i) {
        vtkIdType cellId = psi.getCellId();
        const std::vector<double>& xia = psi.getBegCellParamCoord();
        const std::vector<double>& xib = psi.getEndCellParamCoord();
        double ta = psi.getBegLineParamCoord();
        double tb = psi.getEndLineParamCoord();
        double coeff = psi.getCoefficient();
        std::cout << "test2CellsEdge: seg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                   << '\n';
        psi.next();
    }
    assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);

    ptIds->Delete();
    loc->Delete();
    grid->Delete();
    points->Delete();
}


int main(int argc, char** argv) {

    test1Cell();
    test2Cells();
    test2CellsEdge();

    return 0;
}
