#include <vmtCellLocator.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
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

void testUniformLatLonGrid(int nx, int ny, int numCellsPerBucket) {

    // target points
    double point[3];

    int numCells = nx * ny;
    int numPoints = 4 * numCells;
    double dx = 360. / (double) nx;
    double dy = 180. / (double) ny;

    vtkDoubleArray* coords = vtkDoubleArray::New();
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

    vtkPoints* points = vtkPoints::New();
    points->SetData(coords);

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
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

    // save the grid to file
    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    writer->SetFileName("testCellLocator_latlon.vtk");
    writer->SetInputData(grid);
    writer->Update();
    writer->Delete();

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    cloc->SetDataSet(grid);
    cloc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    cloc->BuildLocator();

    cloc->printBuckets();
    
    // check
    double pBegPtr[] = {0.0, -90.0, 0.0};
    double pEndPtr[] = {360.0, 90.0, 0.0};
    Vec3 pBeg(pBegPtr);
    Vec3 pEnd(pEndPtr);
    std::vector< std::pair<vtkIdType, Vec2> > cellIdLambdas;
    double totLambda;

    vtkIdList* cellIds = vtkIdList::New();
    cloc->FindCellsAlongLine(&pBeg[0], &pEnd[0], 1.e-10, cellIds);
    std::cout << "line " << pBeg << " -> " << pEnd << " intercepts cells:\n";
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        std::cout << cellIds->GetId(i) << ' ';
        if ((i + 1) % 10 == 0) std::cout << '\n';
    }
    std::cout << '\n';

    // across the domain
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (auto& cIdLam : cellIdLambdas) {
        vtkIdType cellId = cIdLam.first;
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        std::cout << "... cellId = " << cellId << " lambda = " << lamIn << " -> " << lamOut << '\n';
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny << ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // the other way across the domain
    pBeg[0] = 360.0; pBeg[1] = -90.0;
    pEnd[0] = 0.0; pEnd[1] = 90.0;    
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        vtkIdType cellId = cIdLam.first;
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        std::cout << "... cellId = " << cellId << " lambda = " << lamIn << " -> " << lamOut << '\n';
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny << ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // along y edge of domain
    pBeg[0] = 0.0; pBeg[1] = -90.0;
    pEnd[0] = 0.0; pEnd[1] = 90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        vtkIdType cellId = cIdLam.first;
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // along x edge of domain
    pBeg[0] = 0.0; pBeg[1] = 90.0;
    pEnd[0] = 360.0; pEnd[1] = 90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // both points are just inside
    pBeg[0] = 0.001; pBeg[1] = -90.0;
    pEnd[0] = 0.001; pEnd[1] = +90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // both points are just outside
    pBeg[0] = -0.001; pBeg[1] = -90.0;
    pEnd[0] = -0.001; pEnd[1] = +90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd, 360.0);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // start point is outside, end point is inside
    pBeg[0] = -0.001; pBeg[1] = -90.0;
    pEnd[0] = +0.001; pEnd[1] = +90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);

    // start point is inside, end point is outside
    pBeg[0] = -0.001; pBeg[1] = -90.0;
    pEnd[0] = +0.001; pEnd[1] = +90.0;
    cellIdLambdas = cloc->findIntersectionsWithLine(pBeg, pEnd);
    totLambda = 0.0;
    for (const auto& cIdLam : cellIdLambdas) {
        double lamIn = cIdLam.second[0];
        double lamOut = cIdLam.second[cIdLam.second.size() - 1];
        totLambda += lamOut - lamIn;
    }
    std::cout << "testUniformLatLonGrid(" << nx << ',' << ny <<  ',' <<
                 numCellsPerBucket << "): pBeg = " << pBeg << " pEnd = " << pEnd <<
                " totLambda = " << totLambda << '\n';
    assert(std::abs(totLambda - 1.0) < 1.e-10);


    // clean up
    cloc->Delete();
    grid->Delete();
    points->Delete();
    coords->Delete();

}


int main(int argc, char** argv) {

    test1Quad(1);
    test1Quad(10);
    test1Quad(1000);
    testUniformLatLonGrid(10, 5, 1);
    testUniformLatLonGrid(10, 5, 2);
    testUniformLatLonGrid(10, 5, 5);
    testUniformLatLonGrid(10, 5, 10);
    testUniformLatLonGrid(10, 5, 20);
    testUniformLatLonGrid(10, 5, 50);

    return 0;
}
