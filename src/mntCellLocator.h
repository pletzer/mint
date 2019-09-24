#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <mntGrid.h>

#ifndef MNT_CELL_LOCATOR
#define MNT_CELL_LOCATOR

/**
 * A class to quickly find a cell in an unstructured grid 
 */

struct CellLocator_t {
    double weights[8]; /* big enough to accommodate quads and hexahedra */
    Grid_t* gridt;
    vtkCellLocator* loc;
    vtkGenericCell* cell;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_new(CellLocator_t** self);

/**
 * Load from file
 * @param fort_filename VTK file name without '\0'
 * @param number of characters in fort_filename
 * @return error code (0 is OK)
 * @note you must construct the object before using this method
 */
extern "C"
int mnt_celllocator_load(CellLocator_t** self, const char* fort_filename, size_t n);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_del(CellLocator_t** self);

/**
 * Set the points (vertices) array pointer
 * @param nVertsPerCell numnber of vertices per cell
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 is OK)
 * @note the caller is responsible for managing the memory of the points array, which
 *       is expected to exist until the grid object is destroyed.
 */
extern "C"
int mnt_celllocator_setPointsPtr(CellLocator_t** self, int nVertsPerCell, size_t ncells, const double points[]);

/**
 * Check if cell areas/volumes are positive
 * @param tol tolerance, area volume must be >= to be valid
 * @param numBadCells number of bad cells (output)
 * @paran badCellids pointer to the array of bad cellIds
 * @return error code (0 is OK)
 * @note pointer badCellIds is invalid if numBadCells is zero
 */
extern "C"
int mnt_celllocator_checkGrid(CellLocator_t** self, double tol, int* numBadCells);

/**
 * Run grid diagnostics
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_runGridDiagnostics(CellLocator_t** self);

/**
 * Build the regridder
 * @param num_cells_per_bucket number of cells per tree node (bucket)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_build(CellLocator_t** self, int num_cells_per_bucket);

/**
 * Save the grid to a VTK file
 * @param filename fortran file name
 * @param n number of characters of fortr_filename
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_dumpGrid(CellLocator_t** self, const char* fort_filename, size_t n);

/**
 * Find the cell and the parametric coordinates
 * @param point target point 
 * @param cellId cell Id (< 0 if not found)
 * @param pcoords parametric coordinates in the unit cell (filled in if found)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], long long* cellId, double pcoords[]);

/**
 * Interpolate the position
 * @param cellId cell Id (input)
 * @param pcoords parametric coordinates in the unit cell (input)
 * @param point target point (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_interpPoint(CellLocator_t** self, long long cellId, const double pcoords[], double point[]);

extern "C"
void mnt_celllocator_printAddress(void* something);

#endif // MNT_CELL_LOCATOR
