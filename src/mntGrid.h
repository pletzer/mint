#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGenericCell.h>

#ifndef MNT_GRID
#define MNT_GRID

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
int mnt_grid_new(Grid_t** self);

/**
 * Destructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_del(Grid_t** self);

/**
 * Set points
 * @param self instance of Grid_t
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_setpoints(Grid_t** self, int nVertsPerCell, int ncells, const double points[]);

/**
 * Get the VTK unstructured grid
 * @param self instance of Grid_t
 * @param grid_ptr prointer to the VTK grid
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr);

/**
 * Load grid from file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory
 */
extern "C"
int mnt_grid_load(Grid_t** self, const char* filename);

/**
 * Load grid to file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename);

#endif // MNT_GRID
