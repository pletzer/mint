#include <vmtCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#include <cassert>
#undef NDEBUG // turn on asserts

void test1Quad(int numCellsPerBucket) {

    double point[3];

    // grid nodes
    vtkDoubleArray* coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    point[0] = 0.0; point[1] = 0.0;; point[2] = 0.0;
    coords->InsertNextTuple(point);
    point[0] = 1.0; point[1] = 0.0;; point[2] = 0.0;
    coords->InsertNextTuple(point);
    point[0] = 1.0; point[1] = 1.0;; point[2] = 0.0;
    coords->InsertNextTuple(point);
    point[0] = 0.0; point[1] = 1.0;; point[2] = 0.0;
    coords->InsertNextTuple(point);

	vtkPoints* points = vtkPoints::New();
    points->SetData(coords);

    // grid
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    grid->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (size_t i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    grid->InsertNextCell(VTK_QUAD, ptIds);
    ptIds->Delete();

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    cloc->SetDataSet(grid);
    cloc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    cloc->BuildLocator();

    // check FindCell

    vtkIdType cellId = -1;
    double tol2 = 1.e-10;
    vtkGenericCell* cell = NULL;
    double pcoords[3];
    double weights[8];

    // point is inside
    point[0] = 0.1; point[1] = 0.2; point[2] = 0.0;
    cellId = cloc->FindCell(point, tol2, cell, pcoords, weights);
    std::cout << "cellId = " << cellId << '\n';
    assert(cellId == 0);
    std::cout << "param coords for point " << point[0] << ',' << point[1] << " are " << pcoords[0] << ',' << pcoords[1] << '\n';
    assert(std::abs(pcoords[0] - 0.1) < tol2);
    assert(std::abs(pcoords[1] - 0.2) < tol2);

    // point is outside
    point[0] = -0.1; point[1] = 0.2; point[2] = 0.0;
    cellId = cloc->FindCell(point, tol2, cell, pcoords, weights);
    assert(cellId < 0);

    // check FindCellsAlongLine
    vtkIdList* cellIds = vtkIdList::New();
    double p0[3], p1[3];

    // start/end points inside the cell
    p0[0] = 0.1; p0[1] = 0.2; p0[2] = 0.0;
    p1[0] = 0.2; p1[1] = 0.3; p1[2] = 0.0;    
    cloc->FindCellsAlongLine(p0, p1, tol2, cellIds);
    assert(cellIds->GetNumberOfIds() == 1);
    assert(cellIds->GetId(0) == 0);

    // no intersection, but ok to detect the cell
    p0[0] =-0.1; p0[1] = 0.0; p0[2] = 0.0;
    p1[0] =-0.1; p1[1] = 1.0; p1[2] = 0.0;    
    cloc->FindCellsAlongLine(p0, p1, tol2, cellIds);
    assert(cellIds->GetNumberOfIds() <= 1);

    cellIds->Delete();

    // clean up
    cloc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();

}

void test1Hex() {
}


int main(int argc, char** argv) {

    test1Quad(1);
    test1Quad(10);
    test1Quad(1000);
    test1Hex();

    return 0;
}
