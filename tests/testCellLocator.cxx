#include <vmtCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
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
    ptIds->InsertNextId(0);
    grid->InsertNextCell(VTK_QUAD, ptIds);
    ptIds->Delete();

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    cloc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    cloc->BuildLocator();

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
