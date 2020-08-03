#include <vector>
#include <string>

#ifndef MNT_NC_DIMENSIONS
#define MNT_NC_DIMENSIONS

/**
 * A class that extracts dimensions from a netcdf variable
 */

struct NcDimensions_t {

  // sizes 
  std::vector<size_t> dims;

};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_new(NcDimensions_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_del(NcDimensions_t** self);

/**
 * Read the dimensions of a variable
 * @param ncid netcdf id
 * @param varid variable id
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_read(NcDimensions_t** self, int ncid, int varid);

/**
 * Get the number of dimensions
 * @param ndims number (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_getNumDims(NcDimensions_t** self, int* ndims);

/**
 * Get a dimension
 * @param len size (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_get(NcDimensions_t** self, int i, size_t* len);

/**
 * Print the attributes
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncdimensions_print(NcDimensions_t** self);


#endif // MNT_NC_DIMENSIONS
