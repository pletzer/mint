#include <vector>
#include <map>
#include <string>

#ifndef MNT_NC_FIELD_WRITE
#define MNT_NC_FIELD_WRITE

/**
 * A class to write and write edge fields to/from netcdf files
 */

struct NcFieldwrite_t {

  // netcdf file handle
  int ncid;

  // variable id
  int varid;

  // variable name
  std::string varName;

  // dimension names and values
  std::vector<std::string> dimNames;
  std::vector<size_t> dimSizes;

  // attributes
  std::map<std::string, std::string> attStr;
  std::map<std::string, int> attInt;
  std::map<std::string, double> attDbl;

  bool defined;
};

/**
 * Constructor
 * @param fileName file name (does not require termination character)
 * @param fileNameLen length of filename string (excluding '\0' if present)
 * @param varName variable name (does not require termination character)
 * @param varNameLen length of varName string (excluding '\0' if present)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_new(NcFieldwrite_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_del(NcFieldwrite_t** self);


/**
 * Net the netcdf variable's number of dimensions
 * @param ndims number (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_setNumDims(NcFieldwrite_t** self, int ndims);


/**
 * Set a netcdf variable's dimension
 * @param dimName name
 * @param dimNameLen number of characters in dimName (excluding '\0')
 * @param dim size
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_setDim(NcFieldwrite_t** self, int iAxis, 
                            const char* dimName, int dimNameLen, size_t dim);

/**
 * Set attribute of type string
 * @param attName attribute name as a flat array
 * @param attNameLen length of attribute name
 * @param attVal attribute value
 * @param attValLen length of attribute value
 */
extern "C"
int mnt_ncfieldwrite_setAttStr(NcFieldwrite_t** self,
                               const char* attName, int attNameLen,
                               const char* attVal, int attValLen);

/**
 * Set attribute of type int
 * @param attName attribute name as a flat array
 * @param attNameLen length of attribute name
 * @param attVal attribute value
 */
extern "C"
int mnt_ncfieldwrite_setAttInt(NcFieldwrite_t** self,
                               const char* attName, int attNameLen,
                               int attVal);


/**
 * Set attribute of type double
 * @param attName attribute name as a flat array
 * @param attNameLen length of attribute name
 * @param attVal attribute value
 */
extern "C"
int mnt_ncfieldwrite_setAttsDbl(NcFieldwrite_t** self,
                                const char* attName, int attNameLen,
                                double attVal);


/**
 * Write the netcdf variable 
 * @param data (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_data(NcFieldwrite_t** self, 
                          const double data[]);

/**
 * Write a slice of the netcdf variable
 * @param startInds0 starting indices for each axis/dimension (0 based)
 * @param counts number of data values to write for each axis/dimension
 * @param data (output)
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldwrite_dataSlice(NcFieldwrite_t** self, 
                               const size_t startInds0[], 
                               const size_t counts[], 
                               const double data[]);



extern "C"
int mnt_ncfieldwrite_define(NcFieldwrite_t** self);


#endif // MNT_NC_FIELD_WRITE
