#include <vector>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct mntRegridEdges_t {
    vtkUnstructuredGrid* srcGrid;
    vtkUnstructuredGrid* dstGrid;
    vtkCellLocator* srcLoc;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_new(mntLatLon_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(mntLatLon_t** self);

/**
 * Set the source grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setSrcGrid(mntLatLon_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the dstination grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setDstGrid(mntLatLon_t** self, vtkUnstructuredGrid* grid);


/**
 * Build the regridder
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_build(mntLatLon_t** self);

/**
 * Load the weights from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_load(mntLatLon_t** self, const std::string& filename);

/**
 * Dump the weights to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dump(mntLatLon_t** self, const std::string& filename);


