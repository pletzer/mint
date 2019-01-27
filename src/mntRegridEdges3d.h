#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <mntGrid.h>

#ifndef MNT_REGRID_EDGES_3D
#define MNT_REGRID_EDGES_3D

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct RegridEdges3d_t {
    vtkUnstructuredGrid* srcGrid;
    vtkUnstructuredGrid* dstGrid;
    vtkCellLocator* srcLoc;
    std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> > weights;
    size_t numSrcCells;
    size_t numDstCells;
    size_t numPointsPerCell;
    size_t numEdgesPerCell;
    Grid_t* srcGridObj;
    Grid_t* dstGridObj;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_new(RegridEdges3d_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_del(RegridEdges3d_t** self);

/**
 * Set the source grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_setSrcGrid(RegridEdges3d_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the destination grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_setDstGrid(RegridEdges3d_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the pointer to the source grid points
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param verts cell vertices as flat array [x0, y0, z0, x1, y1, z1, ...]
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_setSrcPointsPtr(RegridEdges3d_t** self, 
                                    size_t nVertsPerCell, size_t ncells, 
                                    const double verts[]);


/**
 * Set the pointer to the destination grid points
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param verts cell vertices as flat array [x0, y0, z0, x1, y1, z1, ...]
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_setDstPointsPtr(RegridEdges3d_t** self, 
                                    size_t nVertsPerCell, size_t ncells, 
                                    const double verts[]);

/**
 * Build the regridder
 * @param numCellsPerBucket average number of cells per bucket
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_build(RegridEdges3d_t** self, int numCellsPerBucket);

/**
 * Get number of source grid cells
 * @param n number of cells
 */
extern "C"
int mnt_regridedges3d_getNumSrcCells(RegridEdges3d_t** self, int* n);

/**
 * Get number of destination grid cells
 * @param n number of cells
 */
extern "C"
int mnt_regridedges3d_getNumDstCells(RegridEdges3d_t** self, int* n);

/**
 * Get number of edges per cell
 * @param n number
 */
extern "C"
int mnt_regridedges3d_getNumEdgesPerCell(RegridEdges3d_t** self, int* n);

/**
 * Apply interpolation weights to field
 * @param src_data edge centred data on the source grid 
 * @param dst_data edge centred data on the destination grid
 * @return error code (0 is OK)
 * @note edges go anticlockwise
 */
extern "C"
int mnt_regridedges3d_applyWeights(RegridEdges3d_t** self, const double src_data[], double dst_data[]);

/**
 * Load the weights from file
 * @param filename file name
 * @return error code (0 is OK)
 * @note this does not create the object, user must call mnt_regridedges_new prior to this call
 */
extern "C"
int mnt_regridedges3d_load(RegridEdges3d_t** self, const char* filename);

/**
 * Dump the weights to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges3d_dump(RegridEdges3d_t** self, const char* filename);

#endif // MNT_REGRID_EDGES_3D
