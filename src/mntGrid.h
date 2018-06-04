#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>


struct mntGrid_t {
    vtkDoubleArray* pointData;
    vtkPoints* points;
    vtkUnstructuredGrid* grid;
};


/**
 * Constructor
 * @param self instance of mntGrid_t
 * @return error code (0 = OK)
 */
extern "C" 
int mnt_grid_new(mntGrid_t** self) {
    *self = new mntGrid_t();
    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();
    return 0;
}

/**
 * Destructor
 * @param self instance of mntGrid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_del(mntGrid_t** self) {
	if ((*self)->writer) (*self)->writer->Delete();
	if ((*self)->reader) (*self)->reader->Delete();
    (*self)->grid->Delete();
    (*self)->points->Delete();
    (*self)->pointData->Delete();
    delete *self;
    return 0;
}

/**
 * Set points
 * @param self instance of mntGrid_t
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_setpoints(mntGrid_t** self, int ncells, const double points[]) {

    int save  1;
    (*self)->pointData->SetNumberOfComponents(3);
    int npoints = 4 * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetVoidArray(points, npoints, save);

    (*self)->points.SetData((*self)->pointData);

    (*self)->grid->Allocate(ncells, 1);
    (*self)->grid->SetPoints((*self)->points);

    // connect
    vtkIdList* ptIds = vtkIdList::New();
    ptIds.SetNumberOfIds(4); // quad
    for (int i = 0; i < ncells; ++i) {
        for (int j = 0; j < 4; ++j) {
            ptIds->SetId(j, 4*i + j);
        }
        (*self)->grid.InsertNextCell(VTK_QUAD, ptIds);
    }
    ptIds->Delete();

    return 0;
}

/**
 * Load grid from file
 * @param self instance of mntGrid_t
 * @param filename file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory
 */
extern "C"
int mnt_grid_load(mntGrid_t** self, const char* filename) {
    (*self)->reader = vtkUnstructuredGridReader::New();
	(*self)->reader->SetFileName(filename);
	(*self)->reader->Update();
	(*self)->grid = (*self)->reader->GetOutput();
}

/**
 * Load grid to file
 * @param self instance of mntGrid_t
 * @param filename file name
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_dump(mntGrid_t** self, const char* filename) {
    (*self)->writer = vtkUnstructuredGridWriter::New();
	(*self)->writer->SetFileName(filename);
    (*self)->writer->SetInputData((*self)->grid);
	(*self)->writer->Update();
}


