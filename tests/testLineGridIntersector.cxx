#include <mntLineGridIntersector.h>
#include "latLonGrid.h"
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>


void test1Cell() {

    std::cout << "=======  test1Cell ===========\n";

    // a one cell grid, built from scratch
    double v0[] = {0., 0., 0.};
    double v1[] = {1., 0., 0.};
    double v2[] = {1., 1., 0.};
    double v3[] = {0., 1., 0.};
    double v4[] = {0., 0., 1.};
    double v5[] = {1., 0., 1.};
    double v6[] = {1., 1., 1.};
    double v7[] = {0., 1., 1.};
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
    grid->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(8);
    for (vtkIdType i = 0; i < 8; ++i) {
        ptIds->SetId(i, i);
    }
    grid->InsertNextCell(VTK_HEXAHEDRON, ptIds);
    

    LineGridIntersector intersector(grid);

    {
        // line is outiside the cell
        const double p0[] = {1.2, 0.3, 0.};
        const double p1[] = {1.4, 0.3, 0.};
        std::cout << " p0 = "; for (size_t i = 0; i < 3; ++i) std::cout << p0[i] << ',';
        std::cout << " p1 = "; for (size_t i = 0; i < 3; ++i) std::cout << p1[i] << ',';
        std::cout << '\n';
        intersector.setLine(p0, p1);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tintersect t value =  " << tValues[i] << '\n';
        }
        assert(tValues.size() == 0);
    }

    {
        // start/end points are the same
        const double p0[] = {0.2, 0.3, 0.};
        const double p1[] = {0.2, 0.3, 0.};
        std::cout << " p0 = "; for (size_t i = 0; i < 3; ++i) std::cout << p0[i] << ',';
        std::cout << " p1 = "; for (size_t i = 0; i < 3; ++i) std::cout << p1[i] << ',';
        std::cout << '\n';
        intersector.setLine(p0, p1);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tintersect t value =  " << tValues[i] << '\n';
        }
        assert(tValues.size() == 1);
    }


    {
        const double p0[] = {0.2, 0.3, 0.};
        const double p1[] = {0.4, 0.5, 0.};
        std::cout << " p0 = "; for (size_t i = 0; i < 3; ++i) std::cout << p0[i] << ',';
        std::cout << " p1 = "; for (size_t i = 0; i < 3; ++i) std::cout << p1[i] << ',';
        std::cout << '\n';
        intersector.setLine(p0, p1);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tintersect t value =  " << tValues[i] << '\n';
        }
        assert(tValues.size() == 2);
    }

    // line starts on a face
    {
        const double p0[] = {0., 0.5, 0.};
        const double p1[] = {0.4, 0.5, 0.};
        std::cout << " p0 = "; for (size_t i = 0; i < 3; ++i) std::cout << p0[i] << ',';
        std::cout << " p1 = "; for (size_t i = 0; i < 3; ++i) std::cout << p1[i] << ',';
        std::cout << '\n';
        intersector.setLine(p0, p1);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tintersect t value =  " << tValues[i] << '\n';
        }
        assert(tValues.size() == 2);
    }

    {
        const double p0[] = {0., 0.5, 0.};
        const double p1[] = {1.0, 0.5, 0.};
        std::cout << " p0 = "; for (size_t i = 0; i < 3; ++i) std::cout << p0[i] << ',';
        std::cout << " p1 = "; for (size_t i = 0; i < 3; ++i) std::cout << p1[i] << ',';
        std::cout << '\n';
        intersector.setLine(p0, p1);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tintersect t value =  " << tValues[i] << '\n';
        }
        assert(tValues.size() == 2);
    }

    grid->Delete();
    points->Delete();
}


void testLatLon(size_t nElv, size_t nLat, size_t nLon) {

    std::cout << "=======  testLatLon " << nElv << "x" << nLat << "x" << nLon << " ===========\n";

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();

    latLonGrid(nElv, nLat, nLon, grid, points);

    //
    // line segment, elv, lat, lon is fully inside the domain
    // 
    size_t nCases = 11;
    const double paLatLonFullyInside[] = {0.5, 0., -180.,
                               0., 0., -180.,
                               0., -90., -180.,
                               0., -90., -180.,
                               1., 90., 180.,
                               0.5, 90., 180., 
                               0., 90, 180.,
                               1., 0., 180., 
                               1., -90., 180.,
                               1., 90., -180., 
                               1., 90., 0.};
    const double pbLatLonFullyInside[] = {0.5, 0., 180.,
                               0., 0., 180.,
                               0., 90., 180.,
                               1., 90., 180.,
                               0., -90., -180.,
                               0.5, -90., -180.,
                               1., -90., -180., 
                               0., 0., -180.,
                               0., -90, 180., 
                               0., -90., 0.,
                               1., 90., 180.};


    LineGridIntersector intersector(grid);

    for (size_t iCase = 0; iCase < nCases; ++iCase) {

        std::cout << "... segment fully inside: case " << iCase << '\n';
        intersector.setLine(&paLatLonFullyInside[3*iCase], 
                            &pbLatLonFullyInside[3*iCase]);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tx t value =  " << tValues[i] << '\n';
        }
    }

    // line segment is fully outside
    nCases = 11;
    const double paLatLonFullyOutside[] = {2.0, 0., -180.,
                               2., 0., -180.,
                               2., -90., -180.,
                               2., -90., -180.,
                               2., 90., 180.,
                               2.5, 90., 180., 
                               2., 90, 180.,
                               2., 0., 180., 
                               2., -90., 180.,
                               2., 90., -180., 
                               2., 90., 0.};
    const double pbLatLonFullyOutside[] = {2.0, 0., 180.,
                               2., 0., 180.,
                               2., 90., 180.,
                               2., 90., 180.,
                               2., -90., -180.,
                               2.5, -90., -180.,
                               2., -90., -180., 
                               2., 0., -180.,
                               2., -90, 180., 
                               2., -90., 0.,
                               2., 90., 180.};

    for (size_t iCase = 0; iCase < nCases; ++iCase) {
        std::cout << "... segment fully outside: case " << iCase << '\n';
        intersector.setLine(&paLatLonFullyOutside[3*iCase], 
                            &pbLatLonFullyOutside[3*iCase]);
        const std::vector<double>& tValues = intersector.getIntersectionLineParamCoords();
        for (size_t i = 0; i < tValues.size(); ++i) {
            std::cout << "\tx t value =  " << tValues[i] << '\n';
        }
    }


    // clean up
    grid->Delete();
    points->Delete();
}


int main(int argc, char** argv) {

    testLatLon(10, 11, 12);
    testLatLon(1, 11, 12);
    testLatLon(1, 5, 10);
    testLatLon(1, 4, 8);
    testLatLon(2, 4, 4);
    testLatLon(1, 2, 4);
    testLatLon(1, 1, 1);
    test1Cell();

    return 0;
}
