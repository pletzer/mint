#include <mntCellLocator.h>
#include <limits>
#include <cstring>

extern "C"
int mnt_celllocator_new(CellLocator_t** self) {
    *self = new CellLocator_t();
    mnt_grid_new(&(*self)->grid);
    (*self)->loc = vtkCellLocator::New();
    (*self)->cell = vtkGenericCell::New();
    return 0;
}

extern "C"
int mnt_celllocator_load(CellLocator_t** self, const char* fort_filename, size_t n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    char* filename = new char[n + 1];
    strncpy(filename, fort_filename, n);
    filename[n] = '\0';
    // load the grid
    int ier = mnt_grid_load(&(*self)->grid, filename);
    if (ier != 0) {
        return ier;
    }
    return 0;
}

extern "C"
int mnt_celllocator_del(CellLocator_t** self) {
    (*self)->cell->Delete();
    (*self)->loc->Delete();
    mnt_grid_del(&(*self)->grid);
    delete *self;
    return 0;
}

extern "C"
int mnt_celllocator_setpoints(CellLocator_t** self, int nVertsPerCell, vtkIdType ncells, 
		              const double points[]) {
    int ier = mnt_grid_setpoints(&(*self)->grid, nVertsPerCell, ncells, points);
    return ier;
}

extern "C"
int mnt_celllocator_build(CellLocator_t** self, int num_cells_per_bucket) {

    // get the grid
    vtkUnstructuredGrid* grd;
    mnt_grid_get(&(*self)->grid, &grd);

    // build the locator
    (*self)->loc->SetNumberOfCellsPerBucket(num_cells_per_bucket);
    (*self)->loc->SetDataSet(grd);
    (*self)->loc->BuildLocator();

    return 0;
}

extern "C"
int mnt_celllocator_dumpgrid(CellLocator_t** self, const char* filename) {
    int ier = mnt_grid_dump(&(*self)->grid, filename);
    return ier;
}

extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], vtkIdType* cellId, 
		         double pcoords[]) {
    // may need to be adjusted
    const double tol2 = 10 * std::numeric_limits<double>::epsilon(); 
    *cellId = (*self)->loc->FindCell((double*) point, tol2, (*self)->cell, pcoords, (*self)->weights);
    return 0;
}

extern "C"
int mnt_celllocator_interp_point(CellLocator_t** self, vtkIdType cellId, const double pcoords[], double point[]) {
    int subId = 0;
    (*self)->grid->grid->GetCell(cellId)->EvaluateLocation(subId, (double*) pcoords, point, (*self)->weights);
    return 0;
}
