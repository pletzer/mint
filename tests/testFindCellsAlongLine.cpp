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

void test(const std::string& filename, const double p0[], const double p1[], double tol,
            int fixLonAcrossDateline, int averageLonAtPole, int degrees, const char* text) {

    int ier;
    Grid_t* grid = NULL;

    ier = mnt_grid_new(&grid);
    assert(ier == 0);

    // cubed-sphere, radians
    ier = mnt_grid_setFlags(&grid, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);

    ier = mnt_grid_loadFromUgrid2DFile(&grid, filename.c_str());
    assert(ier == 0);

    vtkUnstructuredGrid* ugrid = NULL;
    ier = mnt_grid_get(&grid, &ugrid);
    assert(ier == 0);

    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(ugrid);
    loc->BuildLocator();
    if (degrees == 0) {
        loc->setPeriodicityLengthX(2*M_PI);
    } else {
        // degrees = 1
        loc->setPeriodicityLengthX(360.);
    }

    vtkIdList* cellIds = vtkIdList::New();
    loc->FindCellsAlongLine(p0, p1, tol, cellIds);
    std::cout << text << '\n';
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        std::cout << "cell " << cellIds->GetId(i) << '\n';
    }

    cellIds->Delete();
    loc->Delete();
    ier = mnt_grid_del(&grid);
}


int main(int argc, char** argv) {

    {
        const double p0[] = {270., 0., 0.};
        const double p1[] = {270., 90., 0.};
        test("@CMAKE_SOURCE_DIR@/data/latlon4x2.nc$latlon", p0, p1, 0.1, 0, 0, 1, 
             "latlon4x2 vert edge");
    }
    {
        const double p0[] = {360., -90., 0.};
        const double p1[] = {360., 0., 0.};
        test("@CMAKE_SOURCE_DIR@/data/latlon4x2.nc$latlon", p0, p1, 0.1, 0, 0, 1, 
             "latlon4x2 vert edge at end of domain");
    }
    {
        const double p0[] = {-180., 80., 0.};
        const double p1[] = { 180., 80., 0.};
        test("@CMAKE_SOURCE_DIR@/data/lfric_diag_wind.nc$Mesh2d", p0, p1, 0.1, 1, 1, 1,
            "lfric_diag_wind");
    }
    {
        const double p0[] = {M_PI/2.          , 0., 0.};
        const double p1[] = {M_PI/2. + M_PI/8., 0., 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 0.1, 1, 1, 0,
              "mesh_C4 radians high tol");
        // check effect of tolerance
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 1.e-10, 1, 1, 0,
            "mesh_C4 radians low tol");
    }
    {
        const double p0[] = {1.9634954084936207, 0.36548975596819283, 0.};
        const double p1[] = {1.5707963267948966, 0.39269908169872414, 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", p0, p1, 0.1, 1, 1, 0,
              "mesh_C4 radians 2");
    }

    return 0;
}
