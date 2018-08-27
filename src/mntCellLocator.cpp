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
int mnt_celllocator_setpoints(CellLocator_t** self, int nVertsPerCell, size_t ncells, 
		              const double points[]) {
    int ier = mnt_grid_setpoints(&(*self)->gridt, nVertsPerCell, ncells, points);
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
int mnt_celllocator_rungriddiagnostics(CellLocator_t** self) {
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
int mnt_celllocator_dumpgrid(CellLocator_t** self, const char* fort_filename, size_t n) {
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
int mnt_celllocator_interp_point(CellLocator_t** self, long long cellId, const double pcoords[], double point[]) {
    if (cellId < 0) {
        return 1;
    }
    int subId = 0;
    vtkCell* cell = (*self)->gridt->grid->GetCell(cellId);
    cell->EvaluateLocation(subId, (double*) pcoords, point, (*self)->weights);
    return 0;
}
