#include <mntPolysegmentIter.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellLocator.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>

void test1() {
	// a one cell grid
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
		std::cout << "test1: seg " << i << " cell=" << cellId << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
		                           << " tb=" << tb  << " xib=" << xib[0] << ',' << xib[1] << '\n';
		psi.next();
	}

	ptIds->Delete();
	loc->Delete();
	grid->Delete();
	points->Delete();

}


int main(int argc, char** argv) {

    test1();

    return 0;
}
