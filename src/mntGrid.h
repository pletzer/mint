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

    // vertex raw data
    double* verts;

    // VTK data needed to construct a grid
    vtkDoubleArray* pointData;
    vtkPoints* points;
    vtkUnstructuredGrid* grid;

    vtkUnstructuredGridReader* reader;
    vtkUnstructuredGridWriter* writer;
    std::vector<vtkDoubleArray*> doubleArrays;
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
 * Set the points array pointer
 * @param self instance of Grid_t
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 = OK)
 * @note the caller is responsible for managing the memory of the points array, which
 *       is expected to exist until the grid object is destroyed.
 */
extern "C"
int mnt_grid_setPointsPtr(Grid_t** self, int nVertsPerCell, vtkIdType ncells, const double points[]);

/**
 * Attach data
 * @param self instance of Grid_t
 * @param varname data field name
 * @param nDataPerCell number of components per cell
 * @param points flat array of size ncells*nDataPerCell
 * @return error code (0 = OK)
 * @note this object does own the data, caller is responsible for cleaning the data
 */
extern "C"
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]);

/**
 * Get the VTK unstructured grid
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr);

/**
 * Load grid from a 2D Ugrid file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory
 */
extern "C"
int mnt_grid_loadFrom2DUgrid(Grid_t** self, const char* filename);

/**
 * Load grid from a VTK file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory
 */
extern "C"
int mnt_grid_load(Grid_t** self, const char* filename);

/**
 * Dump the grid to a VTK file
 * @param self instance of Grid_t
 * @param filename file name
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename);

/**
 * Print the content of the grid
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_print(Grid_t** self);

#endif // MNT_GRID
