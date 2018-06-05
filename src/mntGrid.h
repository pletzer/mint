#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGenericCell.h>


struct Grid_t {
    vtkDoubleArray* pointData;
    vtkPoints* points;
    vtkUnstructuredGrid* grid;
    vtkUnstructuredGridReader* reader;
    vtkUnstructuredGridWriter* writer;
};


/**
 * Constructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C" 
int mnt_grid_new(Grid_t** self) {
    *self = new Grid_t();
    (*self)->pointData = NULL;
    (*self)->points = NULL;
    (*self)->grid = NULL;
    (*self)->reader = NULL;
    (*self)->writer = NULL;
    return 0;
}

/**
 * Destructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_del(Grid_t** self) {
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

/**
 * Set points
 * @param self instance of Grid_t
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_setpoints(Grid_t** self, int ncells, const double points[]) {

    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();

    int save = 1;
    (*self)->pointData->SetNumberOfComponents(3);
    int npoints = 4 * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetVoidArray((double*) points, npoints, save);

    (*self)->points->SetData((*self)->pointData);

    (*self)->grid->Allocate(ncells, 1);
    (*self)->grid->SetPoints((*self)->points);

    // connect
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4); // quad
    for (int i = 0; i < ncells; ++i) {
        for (int j = 0; j < 4; ++j) {
            ptIds->SetId(j, 4*i + j);
        }
        (*self)->grid->InsertNextCell(VTK_QUAD, ptIds);
    }
    ptIds->Delete();

    return 0;
}

/**
 * Get the VTK unstructured grid
 * @param self instance of Grid_t
 * @param grid_ptr prointer to the VTK grid
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr) {
    *grid_ptr = (*self)->grid;
    return 0;
}

/**
 * Load grid from file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory
 */
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

/**
 * Load grid to file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename) {
    (*self)->writer = vtkUnstructuredGridWriter::New();
	(*self)->writer->SetFileName(filename);
    (*self)->writer->SetInputData((*self)->grid);
	(*self)->writer->Update();
    return 0;
}


