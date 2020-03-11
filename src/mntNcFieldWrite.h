#include <vector>
#include <map>
#include <string>

#ifndef MNT_NC_FIELD_WRITE
#define MNT_NC_FIELD_WRITE

/**
 * A class to write and write edge fields to/from netcdf files
 */

struct NcFieldWrite_t {

  // netcdf file handle
  int ncid;

  // variable id
  int varid;

  // variable name
  std::string varName;

  // dimension names and values
  std::vector<std::string> dimNames;
  std::vector<size_t> dimSizes;

  bool defined;
  bool append;
};

/**
 * Constructor
 * @param fileName file name (does not require termination character)
 * @param fileNameLen length of filename string (excluding '\0' if present)
 * @param varName variable name (does not require termination character)
 * @param varNameLen length of varName string (excluding '\0' if present)
 * @param append 1 for append, otherwise a new file will be created
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_new(NcFieldWrite_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen, 
                        int append);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_del(NcFieldWrite_t** self);


/**
 * Net the netcdf variable's number of dimensions
 * @param ndims number (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_setNumDims(NcFieldWrite_t** self, int ndims);


/**
 * Set a netcdf variable's dimension
 * @param iAxis axis index (zero-based)
 * @param dimName name
 * @param dimNameLen number of characters in dimName (excluding '\0')
 * @param dim size
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_setDim(NcFieldWrite_t** self, int iAxis, 
                            const char* dimName, int dimNameLen, size_t dim);

/**
 * Write the netcdf variable 
 * @param data (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_data(NcFieldWrite_t** self, 
                          const double data[]);

/**
 * Write a slice of the netcdf variable
 * @param startInds0 starting indices for each axis/dimension (0 based)
 * @param counts number of data values to write for each axis/dimension
 * @param data (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_dataSlice(NcFieldWrite_t** self, 
                               const size_t startInds0[], 
                               const size_t counts[], 
                               const double data[]);



/* private methods */

extern "C"
int mnt_ncfieldwrite_define(NcFieldWrite_t** self);

extern "C"
int mnt_ncfieldwrite_inquire(NcFieldWrite_t** self);


#endif // MNT_NC_FIELD_WRITE
