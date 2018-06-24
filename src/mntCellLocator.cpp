#include <mntCellLocator.h>

extern "C"
int mnt_celllocator_new(CellLocator_t** self) {
    *self = new CellLocator_t();
    mnt_grid_new(&(*self)->grid);
    (*self)->loc = vtkCellLocator::New();
    (*self)->cell = vtkGenericCell::New();
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
int mnt_celllocator_setpoints(CellLocator_t** self, int nVertsPerCell, int ncells, const double points[]) {
    mnt_grid_setpoints(&(*self)->grid, nVertsPerCell, ncells, points);
    return 0;
}

extern "C"
int mnt_celllocator_build(CellLocator_t** self) {

    // get the grid
    vtkUnstructuredGrid* grd;
    mnt_grid_get(&(*self)->grid, &grd);

    // build the locator
    (*self)->loc->SetDataSet(grd);
    (*self)->loc->BuildLocator();

    return 0;
}

extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], vtkIdType* cellId, double pcoords[]) {
    double weights[8];
    const double eps = 1.e-12; // may need to adjust
    *cellId = (*self)->loc->FindCell((double*) point, eps, (*self)->cell, pcoords, weights);
    return 0;
}