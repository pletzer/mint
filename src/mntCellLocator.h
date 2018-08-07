#include <vector>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <mntGrid.h>

#ifndef MNT_CELL_LOCATOR
#define MNT_CELL_LOCATOR

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct CellLocator_t {
    Grid_t* grid;
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
 * Constructor from file
 * @param fort_filename VTK file name without '\0'
 * @param number of characters in fort_filename
 * @return error code (0 is OK)
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
 * Set the points (vertices)
 * @param nVertsPerCell numnber of vertices per cell
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_setpoints(CellLocator_t** self, int nVertsPerCell, vtkIdType ncells, const double points[]);

/**
 * Build the regridder
 * @param num_cells_per_bucket number of cells per tree node (bucket)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_build(CellLocator_t** self, int num_cells_per_bucket);

/**
 * Save the grid to a VTK file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_dumpgrid(CellLocator_t** self, const char* filename);

/**
 * Find the cell and the parametric coordinates
 * @param point target point 
 * @param cellId cell Id (< 0 if not found)
 * @param pcoords parametric coordinates in the unit cell (filled in if found)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], vtkIdType* cellId, double pcoords[]);

/**
 * Interpolate the position
 * @param cellId cell Id (input)
 * @param pcoords parametric coordinates in the unit cell (input)
 * @param point target point (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_celllocator_interp_point(CellLocator_t** self, vtkIdType cellId, const double pcoords[], double point[]);


#endif // MNT_CELL_LOCATOR
