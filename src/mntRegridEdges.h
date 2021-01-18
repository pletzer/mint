
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
 * Constructs a regridding object for edge centred fields
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
 * Set the source grid flags
 * @param self instance of the regridding object
 * @param fixLonAcrossDateline set this to 1 if a periodicty length (360) should be added/subtracted to nodes in order to make each cell as compact as possible
 * @param averageLonAtPole set this to 1 if longitudes at the poles should take the average value of the cell's longitudes
 * @return error code (0 = OK)
 * @note Longitudes are not uniquely defined at the poles
 */
extern "C"
int mnt_regridedges_setSrcGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole);

/**
 * Set the destination grid flags
 * @param self instance of the regridding object
 * @param fixLonAcrossDateline set this to 1 if a periodicty length (360) should be added/subtracted to nodes in order to make each cell as compact as possible
 * @param averageLonAtPole set this to 1 if longitudes at the poles should take the average value of the cell's longitudes
 * @return error code (0 = OK)
 * @note Longitudes are not uniquely defined at the poles
 */
extern "C"
int mnt_regridedges_setDstGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole);

/** 
 * Dump the source grid to a VTK file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_dumpSrcGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength);

/** 
 * Dump the destination grid to a VTK file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_dumpDstGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength);

/** 
 * Inititalize the slice iterator
 * @param self instance of the regridding object
 * @param src_fort_filename src file name (does not require termination character)
 * @param src_nFilenameLength length of src filename string (excluding '\0' if present)
 * @param dst_fort_filename dst file name (does not require termination character)
 * @param dst_nFilenameLength length of dst filename string (excluding '\0' if present)
 * @param append set this to 1 in order to append the data to an exisiting file, 0 otherwise
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param numSlices number of slices (output)
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 * @note Slicing allows one to iterate over non-horizontal dimensions (e.g. height, time, etc.). The same weights 
 *       can then be applied to each additional dimension of the field.
 */
extern "C"
int mnt_regridedges_initSliceIter(RegridEdges_t** self,
                                  const char* src_fort_filename, int src_nFilenameLength,
                                  const char* dst_fort_filename, int dst_nFilenameLength,
                                  int append,
                                  const char* field_name, int nFieldNameLength, 
                                  size_t* numSlices);


/** 
 * Load a slice of the source field from the 2D UGRID file
 * @param self instance of the regridding object
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 * @note Call this method until the return code is != 0 t(no more slices to read)
 * @note Slicing allows one to iterate over non-horizontal dimensions (e.g. height, time, etc.). The same weights 
 *       can then be applied to each additional dimension of the field.
 */
extern "C"
int mnt_regridedges_loadSrcSlice(RegridEdges_t** self, double data[]);


/** 
 * Dump a slice of the destination field slice to 2D UGRID file
 * @param self instance of the regridding object
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 * @note Slicing allows one to iterate over non-horizontal dimensions (e.g. height, time, etc.). The same weights 
 *       can then be applied to each additional dimension of the field.
 */
extern "C"
int mnt_regridedges_dumpDstSlice(RegridEdges_t** self, double data[]);


/** 
 * Increment the slice iterator
 * @param self instance of the regridding object
 * @return error code (0 is OK)
 * @note Slicing allows one to iterate over non-horizontal dimensions (e.g. height, time, etc.). The same weights 
 *       can then be applied to each additional dimension of the field.
 */
extern "C"
int mnt_regridedges_nextSlice(RegridEdges_t** self);


/** 
 * Load a field from a 2D UGRID file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (output)
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 * @note Slicing allows one to iterate over non-horizontal dimensions (e.g. height, time, etc.). The same weights 
 *       can then be applied to each additional dimension of the field.
 */
extern "C"
int mnt_regridedges_loadEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, double data[]);

/** 
 * Dump a field to a 2D UGRID file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param nFilenameLength length of filename string (excluding '\0' if present)
 * @param field_name name of the field
 * @param nFieldNameLength length of field_name string (excluding '\0' if present)
 * @param ndata number of edges and size of data
 * @param data array of size number of unique edges (input)
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_dumpEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, const double data[]);

/** 
 * Load a source grid from a 2D UGRID file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_loadSrcGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/** 
 * Load a destination grid from a 2D UGRID file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_loadDstGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * Build the regridder and compute the regridding weights
 * @param self instance of the regridding object
 * @param numCellsPerBucket average number of cells per bucket
 * @param periodX periodicity length (set to 0 if non-periodic)
 * @param debug 0=no debug info, 1=print debug info, 2=save bad edges in VTK file
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_build(RegridEdges_t** self, int numCellsPerBucket, double periodX, int debug);

/**
 * Get the number of the source grid cells
 * @param self instance of the regridding object
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumSrcCells(RegridEdges_t** self, size_t* n);

/**
 * Get the number of the destination grid cells
 * @param self instance of the regridding object
 * @param n number of cells
 */
extern "C"
int mnt_regridedges_getNumDstCells(RegridEdges_t** self, size_t* n);

/**
 * Get the number of edges per cell
 * @param self instance of the regridding object
 * @param n number (output)
 */
extern "C"
int mnt_regridedges_getNumEdgesPerCell(RegridEdges_t** self, int* n);

/**
 * Get the number of unique edges in the source grid
 * @param self instance of the regridding object
 * @param n number (output)
 * @note Most cells will share edges with their neighbours
 */
extern "C"
int mnt_regridedges_getNumSrcEdges(RegridEdges_t** self, size_t* n);

/**
 * Get the number of unique edges in the destination grid
 * @param self instance of the regridding object
 * @param n number (output)
 * @note Most cells will share edges with their neighbours
 */
extern "C"
int mnt_regridedges_getNumDstEdges(RegridEdges_t** self, size_t* n);

/**
 * Apply the interpolation weights to an edge field with unique edge Ids
 * @param self instance of the regridding object
 * @param src_data edge centred data on the source grid
 * @param dst_data edge centred data on the destination grid
 * @return error code (0 is OK)
 * @note Edges go anticlockwise
 */
extern "C"
int mnt_regridedges_apply(RegridEdges_t** self, 
                          const double src_data[], double dst_data[]);

/**
 * Load the weights from file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param n length of filename string
 * @return error code (0 is OK)
 * @note This does not create the object, user must call mnt_regridedges_new prior to this call
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_loadWeights(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * Dump the weights to file
 * @param self instance of the regridding object
 * @param fort_filename file name (does not require termination character)
 * @param n length of above filename string
 * @return error code (0 is OK)
 * @note Supplying the length of the filename string allows one to call this function from Fortran
 */
extern "C"
int mnt_regridedges_dumpWeights(RegridEdges_t** self, 
                                const char* fort_filename, int n);

/**
 * Print/display the weights
 * @param self instance of the regridding object
 * @return error code (0 is OK)
 */
extern "C"
int mnt_regridedges_print(RegridEdges_t** self);


#endif // MNT_REGRID_EDGES
