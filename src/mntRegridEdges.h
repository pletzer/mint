#include <vector>
#include <map>
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
    std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> > weights;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_new(mntRegridEdges_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(mntRegridEdges_t** self);

/**
 * Set the source grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setSrcGrid(mntRegridEdges_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the dstination grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setDstGrid(mntRegridEdges_t** self, vtkUnstructuredGrid* grid);


/**
 * Build the regridder
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_build(mntRegridEdges_t** self);

/**
 * Load the weights from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_load(mntRegridEdges_t** self, const std::string& filename);

/**
 * Dump the weights to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dump(mntRegridEdges_t** self, const std::string& filename);


