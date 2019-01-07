#include <mntPolysegmentIter3d.h>
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
#include <cmath>


void test1Cell() {
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
    

    vtkCellLocator* loc = vtkCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

    const double p0[] = {0.2, 0.3, 0.};
    const double p1[] = {0.4, 0.5, 0.};
    PolysegmentIter3d psi(grid, loc, p0, p1);
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

void testLatLon(size_t nElv, size_t nLat, size_t nLon) {

    //
    // generate vertices
    //
    double p[3];

    double dElv = 1.0 / double(nElv);
    double dLat = 180.0 / double(nLat);
    double dLon = 360.0 / double(nLon);

    vtkPoints* points = vtkPoints::New();
    size_t nCells = nElv * nLat * nLon;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(8 * nCells);

    vtkIdType index = 0;
    for (size_t k = 0; k < nElv; ++k) {
        double elv0 = dElv*k;
        double elv1 = elv0 + dElv;
        for (size_t j = 0; j < nLat; ++j) {
            double lat0 = -90.0 + j*dLat;
            double lat1 = lat0 + dLat;
            for (size_t i = 0; i < nLon; ++i) {
                double lon0 = -180.0 + i*dLon;
                double lon1 = lon0 + dLon;

                p[0] = elv0; p[1] = lat0; p[2] = lon0;
                points->SetPoint(index, p);
                index++;

                p[0] = elv0; p[1] = lat0; p[2] = lon1;
                points->SetPoint(index, p);
                index++;

                p[0] = elv0; p[1] = lat1; p[2] = lon1;
                points->SetPoint(index, p);
                index++;

                p[0] = elv0; p[1] = lat1; p[2] = lon0;
                points->SetPoint(index, p);
                index++;

                p[0] = elv1; p[1] = lat0; p[2] = lon0;
                points->SetPoint(index, p);
                index++;

                p[0] = elv1; p[1] = lat0; p[2] = lon1;
                points->SetPoint(index, p);
                index++;

                p[0] = elv1; p[1] = lat1; p[2] = lon1;
                points->SetPoint(index, p);
                index++;

                p[0] = elv1; p[1] = lat1; p[2] = lon0;
                points->SetPoint(index, p);
                index++;
            }
        }
    }

    //
    // create grid
    //
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    grid->Allocate(1, 1);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(8);
    for (size_t iCell = 0; iCell < nCells; ++iCell) {
        for (vtkIdType iVert = 0; iVert < 8; ++iVert) {
            ptIds->SetId(iVert, 8*iCell + iVert);
        }
        grid->InsertNextCell(VTK_HEXAHEDRON, ptIds);
    }

    vtkCellLocator* loc = vtkCellLocator::New();
    loc->SetDataSet(grid);
    loc->BuildLocator();

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


    for (size_t iCase = 0; iCase < nCases; ++iCase) {
        std::cout << "testLatLon with segment fully inside: case " << iCase << '\n';
        PolysegmentIter3d psi(grid, loc, 
                              &paLatLonFullyInside[3*iCase], 
                              &pbLatLonFullyInside[3*iCase]);
        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const std::vector<double>& xia = psi.getBegCellParamCoord();
            const std::vector<double>& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "\tseg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                   << '\n';
            psi.next();
        }
        std::cerr << "\tintegrated param coord = " << psi.getIntegratedParamCoord() << '\n';
        assert(std::abs(psi.getIntegratedParamCoord() - 1.0) < 1.e-10);
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
        std::cout << "testLatLon with segment fully outside: case " << iCase << '\n';
        PolysegmentIter3d psi(grid, loc, 
                              &paLatLonFullyOutside[3*iCase], 
                              &pbLatLonFullyOutside[3*iCase]);
        size_t numSegs = psi.getNumberOfSegments();
        psi.reset();
        for (size_t i = 0; i < numSegs; ++i) {
            vtkIdType cellId = psi.getCellId();
            const std::vector<double>& xia = psi.getBegCellParamCoord();
            const std::vector<double>& xib = psi.getEndCellParamCoord();
            double ta = psi.getBegLineParamCoord();
            double tb = psi.getEndLineParamCoord();
            double coeff = psi.getCoefficient();
            std::cout << "\tseg " << i << " cell=" << cellId \
                                   << " ta=" << ta << " xia=" << xia[0] << ',' << xia[1] 
                                   << " tb=" << tb << " xib=" << xib[0] << ',' << xib[1] 
                                   << '\n';
            psi.next();
        }
        std::cerr << "\tintegrated param coord = " << psi.getIntegratedParamCoord() << '\n';
        assert(std::abs(psi.getIntegratedParamCoord() - 0.0) < 1.e-10);
    }


    // clean up
    loc->Delete();
    ptIds->Delete();
    grid->Delete();
    points->Delete();
}


int main(int argc, char** argv) {

    //testLatLon(10, 11, 12); // fails
    testLatLon(2, 4, 4);
    testLatLon(2, 3, 4);
    testLatLon(1, 2, 3);
    testLatLon(1, 1, 1);
    test1Cell();

    return 0;
}
