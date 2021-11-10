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

void test(const std::string& filename, size_t npoints, const double points[]) {

    int ier;
    Grid_t* grid = NULL;

    ier = mnt_grid_new(&grid);
    assert(ier == 0);

    // fix for cubed-sphere, radians
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
    // loc->setPeriodicityLengthX(2*M_PI);
    assert(npoints > 1);

    for (size_t ipoint = 0; ipoint < npoints - 1; ++ipoint) {
        const double* p0 = &points[3*(ipoint + 0)];
        const double* p1 = &points[3*(ipoint + 1)];
        PolysegmentIter psi(ugrid, loc, p0, p1, 2*M_PI);
        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "test: seg " << i << " cell=" << cellId \
                                      << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                      << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                      << '\n';
            psi.next();
        }
        std::cout << "num segments = " << numSegs << " integrated param coord = " << psi.getIntegratedParamCoord() << '\n';
        if (std::abs(psi.getIntegratedParamCoord() - 1.0) >= 1.e-10) {
            char outname[] = "polysegmentIter2.vtk";
            std::cout << "saving grid in file " << outname << '\n';
            ier = mnt_grid_dump(&grid, outname);
        }
        assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);
    }


    loc->Delete();
    ier = mnt_grid_del(&grid);
}

void testUM(const std::string& filename, size_t npoints, const double points[]) {

    int ier;
    Grid_t* grid = NULL;

    ier = mnt_grid_new(&grid);
    assert(ier == 0);

    // fix for cubed-sphere, degrees
    ier = mnt_grid_setFlags(&grid, 0, 0, 1);
    assert(ier == 0);

    ier = mnt_grid_loadFromUgrid2D(&grid, filename.c_str());
    assert(ier == 0);

    vtkUnstructuredGrid* ugrid = NULL;
    ier = mnt_grid_get(&grid, &ugrid);
    assert(ier == 0);

    vmtCellLocator* loc = vmtCellLocator::New();
    loc->SetDataSet(ugrid);
    loc->BuildLocator();
    // loc->setPeriodicityLengthX(2*M_PI);
    assert(npoints > 1);

    for (size_t ipoint = 0; ipoint < npoints - 1; ++ipoint) {
        const double* p0 = &points[3*(ipoint + 0)];
        const double* p1 = &points[3*(ipoint + 1)];
        PolysegmentIter psi(ugrid, loc, p0, p1, 2*M_PI);
        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const Vec3& xia = psi.getBegCellParamCoord();
            const Vec3& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "test: seg " << i << " cell=" << cellId \
                                      << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                      << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                      << '\n';
            psi.next();
        }
        std::cout << "num segments = " << numSegs << " integrated param coord = " << psi.getIntegratedParamCoord() << '\n';
        if (std::abs(psi.getIntegratedParamCoord() - 1.0) >= 1.e-10) {
            char outname[] = "polysegmentIter2.vtk";
            std::cout << "saving grid in file " << outname << '\n';
            ier = mnt_grid_dump(&grid, outname);
        }
        assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);
    }


    loc->Delete();
    ier = mnt_grid_del(&grid);
}

int main(int argc, char** argv) {

    const double eps = 1.73654365e-10;

    // the following fails because of the bending near the pole
    // {
    //     const double points[] = {    0., 1.4, 0.,
    //                              2*M_PI, 1.4, 0.};
    //     test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    // }

    {
        const double points[] = {    0., 1., 0.,
                                 2*M_PI, 1., 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {M_PI/2.              , 0.-eps, 0.,
                                 M_PI/2. + 2.*M_PI/16., 0.+eps, 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0., 1.1, 0.,
                                 1., 1.1, 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0., 0.1, 0.,
                                 2., 0.1, 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0., 1.1, 0.,
                                 2., 1.1, 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0., 1., 0.,
                                 0., -1., 0.,
                                 2*M_PI, -1., 0.};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 3, points);
    }
    {
        const double points[] = {3.9269908169872414, -0.6154797086703874, 0.0,
                                 3.5342917352885173, -0.7458526607730738, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {3.5342917352885173, 0.7458526607730737, 0.0,
                                 3.9269908169872414, 0.6154797086703873, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0.7853981633974483, -1.04089353704597, 0.0,
                                 0.3926990816987242, -0.7458526607730738, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {-0.0, -1.1780972450961724, 0.0,
                                 0.7853981633974483, -1.04089353704597, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {0.0, 1.1780972450961724, 0.0,
                                 0.7853981633974483, 1.04089353704597, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }
    {
        const double points[] = {3.9269908169872414, 1.04089353704597, 0.0,
                                 3.141592653589793, 1.1780972450961724, 0.0};
        test("@CMAKE_SOURCE_DIR@/data/mesh_C4.nc$unit_test", 2, points);
    }    

    return 0;
}
