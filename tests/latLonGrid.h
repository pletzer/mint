#undef NDEBUG // turn on asserts
#include <cmath>
#include <limits> // required by vtkUnstructuredGrid
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>

void latLonGrid(size_t nElv, size_t nLat, size_t nLon, vtkUnstructuredGrid* grid, vtkPoints* points) {

    //
    // generate vertices
    //
    double p[3];

    double dElv = 1.0 / double(nElv);
    double dLat = 180.0 / double(nLat);
    double dLon = 360.0 / double(nLon);

    points->Reset();
    grid->Reset();
    vtkIdList* ptIds = vtkIdList::New();

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
    grid->SetPoints(points);
    grid->Allocate(1, 1);
    ptIds->SetNumberOfIds(8);
    for (size_t iCell = 0; iCell < nCells; ++iCell) {
        for (vtkIdType iVert = 0; iVert < 8; ++iVert) {
            ptIds->SetId(iVert, 8*iCell + iVert);
        }
        grid->InsertNextCell(VTK_HEXAHEDRON, ptIds);
    }

    // clean up
    ptIds->Delete();
}
