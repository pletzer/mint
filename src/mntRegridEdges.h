#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include <mntGrid.h>
#include <mntQuadEdgeIter.h>

#ifndef MNT_REGRID_EDGES
#define MNT_REGRID_EDGES

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct RegridEdges_t {

    // pointers to the src/dst unstructured grid
    vtkUnstructuredGrid* srcGrid;
    vtkUnstructuredGrid* dstGrid;

    // cell locator (octree-based) for fast cell search
    vtkCellLocator* srcLoc;

    // interpolation weights and corresponding src/dst grid cell
    // indices and edges indices
    std::vector<vtkIdType> weightDstCellIds;
    std::vector<int> weightDstEdgeIds;

    std::vector<vtkIdType> weightSrcCellIds;
    std::vector<int> weightSrcEdgeIds;

    std::vector<double> weights;

    size_t numSrcCells;
    size_t numDstCells;
    size_t numPointsPerCell;
    size_t numEdgesPerCell;

    Grid_t* srcGridObj;
    Grid_t* dstGridObj;

    QuadEdgeIter edgeConnectivity;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_new(RegridEdges_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(RegridEdges_t** self);


/** 
 * Load source grid from 2D UGRID
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 */
extern "C"
int mnt_regridedges_loadSrc(RegridEdges_t** self, 
                            const char* fort_filename, int n);

/** 
 * Load destination grid from 2D UGRID
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 */
extern "C"
int mnt_regridedges_loadDst(RegridEdges_t** self, 
                            const char* fort_filename, int n);

/**
 * Set the source grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setSrcGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the destination grid
 * @param grid 
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setDstGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid);

/**
 * Set the pointer to the source grid points
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param verts cell vertices as flat array [x0, y0, z0, x1, y1, z1, ...]
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_setSrcPointsPtr(RegridEdges_t** self, 
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
int mnt_regridedges_setDstPointsPtr(RegridEdges_t** self, 
                                    size_t nVertsPerCell, size_t ncells, 
                                    const double verts[]);

/**
 * Build the regridder
 * @param numCellsPerBucket average number of cells per bucket
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_build(RegridEdges_t** self, int numCellsPerBucket);

/**
 * Get number of source grid cells
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumSrcCells(RegridEdges_t** self, int* n);

/**
 * Get number of destination grid cells
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumDstCells(RegridEdges_t** self, int* n);

/**
 * Get number of edges per cell
 * @param n number
 */
extern "C"
int mnt_regridedges_getNumEdgesPerCell(RegridEdges_t** self, int* n);

/**
 * Apply interpolation weights to field
 * @param src_data edge centred data on the source grid 
 * @param dst_data edge centred data on the destination grid
 * @return error code (0 is OK)
 * @note edges go anticlockwise
 */
extern "C"
int mnt_regridedges_applyWeights(RegridEdges_t** self, const double src_data[], double dst_data[]);

/**
 * Load the weights from file
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note this does not create the object, user must call mnt_regridedges_new prior to this call
 */
extern "C"
int mnt_regridedges_load(RegridEdges_t** self, 
                         const char* fort_filename, int n);

/**
 * Dump the weights to file
 * @param fort_filename file name (does not require termination character)
 * @param n length of above filename string
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dump(RegridEdges_t** self, 
                         const char* fort_filename, int n);

/**
 * Print the weights
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_print(RegridEdges_t** self);


#endif // MNT_REGRID_EDGES
