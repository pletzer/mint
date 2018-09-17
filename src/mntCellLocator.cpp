#include <mntCellLocator.h>
#include <limits>
#include <cstring>
#include <string>

extern "C"
int mnt_celllocator_new(CellLocator_t** self) {
    *self = new CellLocator_t();
    mnt_grid_new(&(*self)->gridt);
    (*self)->loc = vtkCellLocator::New();
    (*self)->cell = vtkGenericCell::New();
    return 0;
}

extern "C"
int mnt_celllocator_load(CellLocator_t** self, const char* fort_filename, size_t n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    // load the grid
    int ier = mnt_grid_load(&(*self)->gridt, filename.c_str());
    if (ier != 0) {
        return ier;
    }
    return 0;
}

extern "C"
int mnt_celllocator_del(CellLocator_t** self) {
    (*self)->cell->Delete();
    (*self)->loc->Delete();
    mnt_grid_del(&(*self)->gridt);
    delete *self;
    return 0;
}

extern "C"
int mnt_celllocator_setPointsPtr(CellLocator_t** self, int nVertsPerCell, size_t ncells, 
		                         const double points[]) {
    int ier = mnt_grid_setPointsPtr(&(*self)->gridt, nVertsPerCell, ncells, points);
    return ier;
}

extern "C"
int mnt_celllocator_build(CellLocator_t** self, int num_cells_per_bucket) {

    // get the grid
    vtkUnstructuredGrid* grd;
    mnt_grid_get(&(*self)->gridt, &grd);

    // build the locator
    (*self)->loc->SetNumberOfCellsPerBucket(num_cells_per_bucket);
    (*self)->loc->SetDataSet(grd);
    (*self)->loc->BuildLocator();

    return 0;
}

extern "C"
int mnt_celllocator_checkGrid(CellLocator_t** self, double tol, int* numBadCells) {
    int res = 0;
    vtkIdType ncells = (*self)->gridt->grid->GetNumberOfCells();
    if (ncells == 0) return res;
    vtkIdType npts = (*self)->gridt->grid->GetNumberOfPoints();
    vtkPoints* points = (*self)->gridt->grid->GetPoints();
    vtkCell* cell;
    double vert0[3], vert1[3], vert2[3], vert3[3];
    double a[3], b[3], c[3];
    *numBadCells = 0;
    // check the cell type
    if ((*self)->gridt->grid->GetCell(0)->GetCellType() == VTK_QUAD) {
        for (vtkIdType i = 0; i < ncells; ++i) {
            int numBad = 0;
            vtkIdList* ptIds = (*self)->gridt->grid->GetCell(i)->GetPointIds();
            for (vtkIdType j = 1; j < ptIds->GetNumberOfIds() - 1; j += 2) {
                vtkIdType k0 = ptIds->GetId(j);
                vtkIdType k1 = ptIds->GetId(j + 1);
                vtkIdType k2 = ptIds->GetId(j - 1);
                points->GetPoint(k0, vert0);
                points->GetPoint(k1, vert1);
                points->GetPoint(k2, vert2);
                a[0] = vert1[0] - vert0[0]; a[1] = vert1[1] - vert0[1];
                b[0] = vert2[0] - vert0[0]; b[1] = vert2[1] - vert0[1];
                double area = a[0]*b[1] - a[1]*b[0];
                if (area < std::abs(tol)) {
                    numBad++;
                }
            }
            if (numBad > 0) *numBadCells++;
        }
    }
    else if ((*self)->gridt->grid->GetCell(0)->GetCellType() == VTK_HEXAHEDRON) {
        for (vtkIdType i = 0; i < ncells; ++i) {
            int numBad = 0;
            vtkIdList* ptIds = (*self)->gridt->grid->GetCell(i)->GetPointIds();
            for (vtkIdType j = 1; j < ptIds->GetNumberOfIds() - 1; j += 2) {
                vtkIdType k0 = ptIds->GetId(j);
                vtkIdType k1 = ptIds->GetId(j + 1);
                vtkIdType k2 = ptIds->GetId(j - 1);
                vtkIdType k3 = ptIds->GetId((j + 4) % 8);
                points->GetPoint(k0, vert0);
                points->GetPoint(k1, vert1);
                points->GetPoint(k2, vert2);
                points->GetPoint(k3, vert3);
                a[0] = vert1[0] - vert0[0]; a[1] = vert1[1] - vert0[1]; a[2] = vert1[2] - vert0[2];
                b[0] = vert2[0] - vert0[0]; b[1] = vert2[1] - vert0[1]; b[2] = vert2[2] - vert0[2];
                c[0] = vert3[0] - vert0[0]; c[1] = vert3[1] - vert0[1]; c[2] = vert3[2] - vert0[2];
                double volume = (a[1]*b[2] - a[2]*b[1])*c[0] 
                              + (a[2]*b[0] - a[0]*b[2])*c[1]
                              + (a[0]*b[1] - a[1]*b[0])*c[2];
                if (volume < std::abs(tol)) {
                    numBad++;
                }
            }
            if (numBad > 0) *numBadCells++;
        }
    }
    else {
        res = 1;
    }
    return res;
}

extern "C"
int mnt_celllocator_runGridDiagnostics(CellLocator_t** self) {
    vtkIdType ncells = (*self)->gridt->grid->GetNumberOfCells();
    vtkIdType npts = (*self)->gridt->grid->GetNumberOfPoints();
    std::cout << "No of cells = " << ncells << " no of points = " << npts << '\n';
    vtkPoints* points = (*self)->gridt->grid->GetPoints();
    vtkCell* cell;
    double vert[3];
    for (vtkIdType i = 0; i < ncells; ++i) {
        std::cout << "cell " << i << " verts ";
        cell = (*self)->gridt->grid->GetCell(i);
        vtkIdList* ptIds = cell->GetPointIds();
        for (vtkIdType j = 0; j < ptIds->GetNumberOfIds(); ++j) {
            vtkIdType k = ptIds->GetId(j);
            points->GetPoint(k, vert);
            std::cout << "\t" << k << " (";
            for (int el = 0; el < 3; ++el) {
                std::cout << vert[el] << ',';
            }
            std::cout << ") ";
        }
        std::cout << '\n';
    }
    return 0;
}


extern "C"
int mnt_celllocator_dumpGrid(CellLocator_t** self, const char* fort_filename, size_t n) {
    // copy fortran string into c string
    std::string filename = std::string(fort_filename, n);
    int ier = mnt_grid_dump(&(*self)->gridt, filename.c_str());
    return ier;
}

extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], long long* cellId, 
		                 double pcoords[]) {
    // may need to be adjusted
    const double tol2 = 10 * std::numeric_limits<double>::epsilon(); 
    *cellId = (*self)->loc->FindCell((double*) point, tol2, (*self)->cell, pcoords, (*self)->weights);
    return 0;
}

extern "C"
int mnt_celllocator_interpPoint(CellLocator_t** self, long long cellId, const double pcoords[], double point[]) {
    if (cellId < 0) {
        return 1;
    }
    int subId = 0;
    vtkCell* cell = (*self)->gridt->grid->GetCell(cellId);
    cell->EvaluateLocation(subId, (double*) pcoords, point, (*self)->weights);
    return 0;
}

extern "C"
void mnt_celllocator_printAddress(void* something) {
    printf("address is %lld (0x%16zx)\n", (long long) something, (size_t) something);
}
