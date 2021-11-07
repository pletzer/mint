#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <limits> // required by vtkUnstructuredGrid
#include <mntPolysegmentIter.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <iostream>

void test1(const std::string& filename, const double p0[], const double p1[], double tol) {

    int ier;
    Grid_t* grid = NULL;

    ier = mnt_grid_new(&grid);
    assert(ier == 0);

    // cubed-sphere, radians
    ier = mnt_grid_setFlags(&grid, 1, 1, 0);
    assert(ier == 0);

    ier = mnt_grid_loadFromUgrid2D(&grid, filename.c_str());
    assert(ier == 0);

    vtkUnstructuredGrid* ugrid = NULL;
    ier = mnt_grid_get(&grid, &ugrid);
    assert(ier == 0);

    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(ugrid);
    loc->BuildLocator();
    loc->setPeriodicityLengthX(2*M_PI);

    vtkIdList* cellIds = vtkIdList::New();
    loc->FindCellsAlongLine(p0, p1, tol, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        std::cout << "cell " << cellIds->GetId(i) << '\n';
    }

    cellIds->Delete();
    loc->Delete();
    ier = mnt_grid_del(&grid);
}


int main(int argc, char** argv) {

    {
        const double p0[] = {M_PI/2.          , 0., 0.};
        const double p1[] = {M_PI/2. + M_PI/8., 0., 0.};
        test1("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 0.1);
        test1("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 1.e-10);
    }
    {
        const double p0[] = {1.9634954084936207, 0.36548975596819283, 0.};
        const double p1[] = {1.5707963267948966, 0.39269908169872414, 0.};
        test1("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 1.e-10);
    }

    return 0;
}
