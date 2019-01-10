#undef NDEBUG // turn on asserts
#include <vtkUnstructuredGrid.h>
#include <mvtkCellLocator.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkGenericCell.h>
#include <iostream>
#include <cmath>

void testLatLon(size_t nElv, size_t nLat, size_t nLon, const double pa[], const double pb[]) {

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
    vtkIdType cellId = 0;
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

                std::cout << "cell " << cellId << ": elv = " << elv0 << " -> " << elv1 << 
                             " lat = " << lat0 << " -> " << lat1 << 
                             " lon = " << lon0 << " -> " << lon1 << '\n';

                cellId++;
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
        std::cout << "inserting cell " << iCell;
        for (vtkIdType iVert = 0; iVert < 8; ++iVert) {
          vtkIdType vi = ptIds->GetId(iVert);
          double vrtx[3];
          points->GetPoint(vi, vrtx);
          std::cout << "\t" << iVert << " vert index: " << vi << 
          " vert coords: " << vrtx[0] << ',' << vrtx[1] << ',' << vrtx[2] << '\n';
        }
    }

    mvtkCellLocator* loc = mvtkCellLocator::New();
    loc->SetNumberOfCellsPerBucket(1);
    loc->SetDataSet(grid);
    loc->BuildLocator();

    // find all the cells allong the line
    vtkIdList* cellIds = vtkIdList:: New();
    double tol = 0.001;
    loc->FindCellsAlongLine((double*) pa, (double*) pb, tol, cellIds);
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
      std::cout << "line intersects with cell " << cellIds->GetId(i) << '\n';
    }

    // clean up
    cellIds->Delete();
    loc->Delete();
    ptIds->Delete();
    grid->Delete();
    points->Delete();
}


int main(int argc, char** argv) {

    double pa[] = {0.1, -90.,-180.};
    double pb[] = {0.1, 90., 180.};
    testLatLon(2, 4, 4, pa, pb);

    return 0;
}
