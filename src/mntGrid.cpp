#include <mntGrid.h>
#include <vtkCellData.h>

extern "C" 
int mnt_grid_new(Grid_t** self) {

    *self = new Grid_t();
    (*self)->pointData = NULL;
    (*self)->points = NULL;
    (*self)->grid = NULL;
    (*self)->reader = NULL;
    (*self)->writer = NULL;
    (*self)->doubleArrays.resize(0);
    return 0;
}

extern "C"
int mnt_grid_del(Grid_t** self) {

    for (size_t i = 0; i < (*self)->doubleArrays.size(); ++i) {
        (*self)->doubleArrays[i]->Delete();
    }
    if ((*self)->writer) (*self)->writer->Delete();
    if ((*self)->reader) {
        (*self)->reader->Delete();
    }
    else {
        if ((*self)->grid) (*self)->grid->Delete();
    }
    if ((*self)->points) (*self)->points->Delete();
    if ((*self)->pointData) (*self)->pointData->Delete();
    delete *self;
    return 0;
}

extern "C"
int mnt_grid_setpoints(Grid_t** self, int nVertsPerCell, vtkIdType ncells, const double points[]) {

    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();

    int save = 1;
    int npoints = nVertsPerCell * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetNumberOfComponents(3);
    (*self)->pointData->SetVoidArray((double*) points, npoints*3, save);

    (*self)->points->SetData((*self)->pointData);

    (*self)->grid->Allocate(ncells, 1);
    (*self)->grid->SetPoints((*self)->points);

    int cellType = -1;
    if (nVertsPerCell == 4) {
        cellType = VTK_QUAD;
    }
    else if (nVertsPerCell == 8) {
        cellType = VTK_HEXAHEDRON;
    }
    else {
        // error
        return 1;
    }

    // connect
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(nVertsPerCell);
    for (int i = 0; i < ncells; ++i) {
        for (int j = 0; j < nVertsPerCell; ++j) {
            ptIds->SetId(j, nVertsPerCell*i + j);
        }
        (*self)->grid->InsertNextCell(cellType, ptIds);
    }
    ptIds->Delete();

    return 0;
}

extern "C"
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]) {

    if (!(*self)->grid) {
        return 1;
    }

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();

    vtkDoubleArray* vtkdata = vtkDoubleArray::New();
    vtkdata->SetName(varname);
    vtkdata->SetNumberOfTuples(ncells);
    vtkdata->SetNumberOfComponents(nDataPerCell);
    int save = 1;
    vtkdata->SetVoidArray((double*) data, ncells*nDataPerCell, save);

    // store
    (*self)->doubleArrays.push_back(vtkdata);

    // add to the grid
    (*self)->grid->GetCellData()->AddArray(vtkdata);

    return 0;
}


extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr) {
    *grid_ptr = (*self)->grid;
    return 0;
}

extern "C"
int mnt_grid_load(Grid_t** self, const char* filename) {
    if ((*self)->grid) {
        (*self)->grid->Delete();
    }
    (*self)->reader = vtkUnstructuredGridReader::New();
    (*self)->reader->SetFileName(filename);
    (*self)->reader->Update();
    (*self)->grid = (*self)->reader->GetOutput();
    return 0;
}

extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename) {
    (*self)->writer = vtkUnstructuredGridWriter::New();
    (*self)->writer->SetFileName(filename);
    (*self)->writer->SetInputData((*self)->grid);
    (*self)->writer->Update();
    return 0;
}
