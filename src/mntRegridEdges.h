
#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntQuadEdgeIter.h>
#include <mntNcAttributes.h>
#include <mntNcFieldRead.h>
#include <mntNcFieldWrite.h>
#include <mntMultiArrayIter.h>

#ifndef MNT_REGRID_EDGES
#define MNT_REGRID_EDGES

/**
 * @brief Edge-centred field regridding
 *
 * A class to compute the regridding weights of an edge-centred field
 */

struct RegridEdges_t {

    /** pointer to the source VTK unstructured grid */
    vtkUnstructuredGrid* srcGrid;

    /** pointer to the destination VTK unstructured grid */
    vtkUnstructuredGrid* dstGrid;

    /** cell locator (octree-based) for fast cell search */
    vmtCellLocator* srcLoc;

    /** interpolation weights */
    std::vector<double> weights;

    /** destination cell Ids of the weights */
    std::vector<vtkIdType> weightDstCellIds;

    /** destination face indices of the weights */
    std::vector<int> weightDstFaceEdgeIds;

    /** source cell Ids of the weights */
    std::vector<vtkIdType> weightSrcCellIds;

    /** source face indices of the weights */
    std::vector<int> weightSrcFaceEdgeIds;

    /** number of points per cell */
    size_t numPointsPerCell;

    /** number of edges per cell */
    size_t numEdgesPerCell;

    /** destination grid object */
    Grid_t* dstGridObj;

    /** source grid object */
    Grid_t* srcGridObj;

    /** edge connectivity iterator */
    QuadEdgeIter edgeConnectivity;

    /** number of dimensions */
    int ndims;

    /** netcdf dimnension names */
    std::vector<std::string> dimNames;

    /** netcdf file source data reader */
    NcFieldRead_t* srcReader;

    /** netcdf source file Id */
    int srcNcid;

    /** netcdf variable Id */
    int srcVarid;

    /** netcdf source data dimensions */
    std::vector<size_t> srcDims;

    /** netcdf source data counts */
    std::vector<size_t> srcCounts;

    /** netcdf destination data writer */
    NcFieldWrite_t* dstWriter;

    /** netcdf destination data dimensions */
    std::vector<size_t> dstDims;

    /** netcdf destination data counts */
    std::vector<size_t> dstCounts;

    /** multi-array iterator */
    MultiArrayIter_t* mai;

    /** start indices for multi-array iterator */
    std::vector<size_t> startIndices;

};

/**
 * Constructs a regridding object for edge centred fields. 
 * @param self instance of the regridding object
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_new(RegridEdges_t** self);

/**
 * Destructor
 * @param self this instance
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_del(RegridEdges_t** self);

/**
 * Set source grid flags
 * @param self this instance
 * @param fixLonAcrossDateline set this to 1 if periodicty length should be added/subtracted to nodes in order to make the cell as compact as possible
 * @param averageLonAtPole set this to 1 if longitudes at the poles should take the average value of the node cell node's longitudes
 * @return error code (0 = OK)
 */
extern "C"
int mnt_regridedges_setSrcGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole);

/**
 * Set destination grid flags
 * @param self this instance
 * @param fixLonAcrossDateline set this to 1 if periodicty length should be added/subtracted to nodes in order to make the cell as compact as possible
 * @param averageLonAtPole set this to 1 if longitudes at the poles should take the average value of the node cell node's longitudes
 * @return error code (0 = OK)
 */
extern "C"
int mnt_regridedges_setDstGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole);

/** 
 * Dump source grid to VTK file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dumpSrcGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength);

/** 
 * Dump destination grid to VTK file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dumpDstGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength);

/** 
 * Inititalize source slice iterator
 * @param self this instance
 * @param src_fort_filename src file name (does not require termination character)
 * @param src_nFilenameLength length of src filename string (excluding '\0' if present)
 * @param dst_fort_filename dst file name (does not require termination character)
 * @param dst_nFilenameLength length of dst filename string (excluding '\0' if present)
 * @param append set this to 1 in order to append the data to an exisiting file, 0 otherwise
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param numSlices number of slices (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_initSliceIter(RegridEdges_t** self,
                                  const char* src_fort_filename, int src_nFilenameLength,
                                  const char* dst_fort_filename, int dst_nFilenameLength,
                                  int append,
                                  const char* field_name, int nFieldNameLength, 
                                  size_t* numSlices);


/** 
 * Load a slice of a source field from 2D UGRID file and increment iterator
 * @param self this instance
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 * @note call this method until the return code is != 0 to read each slice
 */
extern "C"
int mnt_regridedges_loadSrcSlice(RegridEdges_t** self, double data[]);


/** 
 * Dump slice of destination field slice to 2D UGRID file
 * @param self this instance
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dumpDstSlice(RegridEdges_t** self, double data[]);


/** 
 * Increment the slice iterator
 * @param self this instance
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_nextSlice(RegridEdges_t** self);


/** 
 * Load field from 2D UGRID file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_loadEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, double data[]);

/** 
 * Dump field to 2D UGRID file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (input)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dumpEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, const double data[]);

/** 
 * Load source grid from 2D UGRID file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_loadSrcGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/** 
 * Load destination grid from 2D UGRID file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_loadDstGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * Build the regridder
 * @param self this instance
 * @param numCellsPerBucket average number of cells per bucket
 * @param periodX periodicity length (set to 0 if non-periodic)
 * @param debug 0=no debug info, 1=print debug info, 2=save bad edges in VTK file
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_build(RegridEdges_t** self, int numCellsPerBucket, double periodX, int debug);

/**
 * Get number of source grid cells
 * @param self this instance
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumSrcCells(RegridEdges_t** self, size_t* n);

/**
 * Get number of destination grid cells
 * @param self this instance
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumDstCells(RegridEdges_t** self, size_t* n);

/**
 * Get number of edges per cell
 * @param self this instance
 * @param n number (output)
 */
extern "C"
int mnt_regridedges_getNumEdgesPerCell(RegridEdges_t** self, int* n);

/**
 * Get number of unique edges in the source grid
 * @param self this instance
 * @param n number (output)
 */
extern "C"
int mnt_regridedges_getNumSrcEdges(RegridEdges_t** self, size_t* n);

/**
 * Get number of unique edges in the destination grid
 * @param self this instance
 * @param n number (output)
 */
extern "C"
int mnt_regridedges_getNumDstEdges(RegridEdges_t** self, size_t* n);

/**
 * Apply interpolation weights to edge field with unique edge Ids
 * @param self this instance
 * @param src_data edge centred data on the source grid
 * @param dst_data edge centred data on the destination grid
 * @return error code (0 is OK)
 * @note edges go anticlockwise
 */
extern "C"
int mnt_regridedges_apply(RegridEdges_t** self, 
                          const double src_data[], double dst_data[]);

/**
 * Load the weights from file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note this does not create the object, user must call mnt_regridedges_new prior to this call
 */
extern "C"
int mnt_regridedges_loadWeights(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * Dump the weights to file
 * @param self this instance
 * @param fort_filename file name (does not require termination character)
 * @param n length of above filename string
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_dumpWeights(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * @param self this instance
 * Print the weights
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_print(RegridEdges_t** self);


#endif // MNT_REGRID_EDGES
