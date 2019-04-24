#include <vector>
#include <map>
#include <string>
#include <mntUgrid2D.h>

#ifndef MNT_REGRID_EDGES_2D
#define MNT_REGRID_EDGES_2D

/**
 * A class to compute the regridding weights of an edge-centred field
 */

struct RegridEdges2D_t {

    // pointers to the src/dst unstructured grid
    Ugrid2D srcGrid;
    Ugrid2D dstGrid;

    // weights
    std::vector<long long> weightSrcEdgeIds;
    std::vector<long long> weightDstEdgeIds;
    std::vector<double> weights;

    size_t numSrcEdges;
    size_t numDstEdges;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges2d_new(RegridEdges2D_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges2d_del(RegridEdges2D_t** self);


/** 
 * Load field from 2D UGRID file
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (output)
 */
extern "C"
int mnt_regridedges2d_loadEdgeField(RegridEdges2D_t** self,
                                    const char* fort_filename, int nFilenameLength,
                                    const char* field_name, int nFieldNameLength,
                                    size_t ndata, double data[]);

/** 
 * Dump field from 2D UGRID file
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (input)
 */
extern "C"
int mnt_regridedges2d_dumpEdgeField(RegridEdges2D_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, const double data[]);

/** 
 * Load source grid from 2D UGRID file
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 */
extern "C"
int mnt_regridedges2d_loadSrcGrid(RegridEdges2D_t** self, 
                                const char* fort_filename, int n);

/** 
 * Load destination grid from 2D UGRID file
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 */
extern "C"
int mnt_regridedges2d_loadDstGrid(RegridEdges2D_t** self, 
                                const char* fort_filename, int n);

/**
 * Get number of unique edges in the source grid
 * @param n number (output)
 */
extern "C"
int mnt_regridedges2d_getNumSrcEdges(RegridEdges2D_t** self, size_t* n);

/**
 * Get number of unique edges in the destination grid
 * @param n number (output)
 */
extern "C"
int mnt_regridedges2d_getNumDstEdges(RegridEdges2D_t** self, size_t* n);


/**
 * Apply interpolation weights to edge field with unique edge Ids
 * @param src_data edge centred data on the source grid
 * @param numDstEdges number of destination grid edges
 * @param dst_data edge centred data on the destination grid
 * @return error code (0 is OK)
 * @note edges go anticlockwise
 */
extern "C"
int mnt_regridedges2d_apply(RegridEdges2D_t** self, 
                            const double src_data[], double dst_data[]);

/**
 * Load the weights from file
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note this does not create the object, user must call mnt_regridedges_new prior to this call
 */
extern "C"
int mnt_regridedges2d_loadWeights(RegridEdges2D_t** self, 
                                const char* fort_filename, int n);

/**
 * Dump the weights to file
 * @param fort_filename file name (does not require termination character)
 * @param n length of above filename string
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges2d_dumpWeights(RegridEdges2D_t** self, 
                                const char* fort_filename, int n);

/**
 * Print the weights
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges2d_print(RegridEdges2D_t** self);


#endif // MNT_REGRID_EDGES_2D
