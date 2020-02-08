#include <vector>
#include <map>
#include <string>

#ifndef MNT_NC_FIELD_READ
#define MNT_NC_FIELD_READ

/**
 * A class to read and write edge fields to/from netcdf files
 */

struct NcFieldRead_t {

  // netcdf file handle
  int ncid;

  // variable id
  int varid;

  // string attributes
  std::map<std::string, std::string> strAttr;

  // dimension names and values
  std::vector<std::string> dimNames;
  std::vector<size_t> dimSizes;
  std::vector<int> dimIds;
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
int mnt_ncfieldread_new(NcFieldRead_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncfieldread_del(NcFieldRead_t** self);


extern "C"
int mnt_ncfieldread_getNumDims(NcFieldRead_t** self, int* ndims);


extern "C"
int mnt_ncfieldread_getDimName(NcFieldRead_t** self, int iAxis, char* dimName, int dimLen);


extern "C"
int mnt_ncfieldread_getDim(NcFieldRead_t** self, int iAxis, size_t* dim);

extern "C"
int mnt_ncfieldread_readData(NcFieldRead_t** self, 
                             double data[]);

extern "C"
int mnt_ncfieldread_readDataSlice(NcFieldRead_t** self, 
                                  const size_t* startInds0, 
                                  const size_t* counts, 
                                  double data[]);



#endif // MNT_NC_FIELD_READ
