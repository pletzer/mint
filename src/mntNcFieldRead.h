#include "mntLIBRARY_API.h"
#include <vector>
#include <map>
#include <string>

#ifndef MNT_NC_FIELD_READ
#define MNT_NC_FIELD_READ


/**
 * Opens a NetCDF file and gets the Id of the file and the variable
 * @param filename NetCDF file name
 * @param varname variable name
 * @param ncid file Id (ouput)
 * @param varid variable Id (ouput)
 * @return 0 = no error, 1 = could not open the file, 2 = could not find the variable
 * @note use mnt_closeNc to close the file
 */
LIBRARY_API
int mnt_openNc(const char* filename, const char* varname, int* ncid, int* varid);

/**
 * Closes a NetCDF file
 * @param ncid file Id
 * @return 0 = no error
 * @note use mnt_openNc to open the file
 */
LIBRARY_API
int mnt_closeNc(int ncid);

/**
 * A class to read and edge fields from netcdf files
 */

struct NcFieldRead_t {

  // netcdf file handle
  int ncid;

  // variable id
  int varid;

  // dimension names and values
  std::vector<std::string> dimNames;
  std::vector<std::size_t> dimSizes;
  std::vector<int> dimIds;
};

/**
 * Constructor
 * @param ncid netcdf file Id
 * @param varid netcdf variable Id
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_new(NcFieldRead_t** self, int ncid, int varid);

/**
 * Destructor
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_del(NcFieldRead_t** self);


/**
 * Get the netcdf variable's number of dimensions
 * @param ndims number (output)
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_getNumDims(NcFieldRead_t** self, int* ndims);


/**
 * Get a netcdf variable's dimension name
 * @param dimName name (output)
 * @param dimNameLen number of characters in dimName (excluding '\0')
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_getDimName(NcFieldRead_t** self, int iAxis, char* dimName, int dimNameLen);


/**
 * Get a netcdf variable dimension's size 
 * @param dim number (output)
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_getDim(NcFieldRead_t** self, int iAxis, std::size_t* dim);

/**
 * Read the netcdf variable 
 * @param data (output)
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_data(NcFieldRead_t** self, 
                         double data[]);

/**
 * Read a slice of the netcdf variable
 * @param startInds0 starting indices for each axis/dimension (0 based)
 * @param counts number of data values to read for each axis/dimension
 * @param data (output)
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_ncfieldread_dataSlice(NcFieldRead_t** self, 
                              const std::size_t* startInds0, 
                              const std::size_t* counts, 
                              double data[]);



#endif // MNT_NC_FIELD_READ
