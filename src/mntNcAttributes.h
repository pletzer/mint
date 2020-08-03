#include <vector>
#include <map>
#include <string>

#ifndef MNT_NC_ATTRIBUTES
#define MNT_NC_ATTRIBUTES

/**
 * A class that stores netcdf attributes
 */

struct NcAttributes_t {

  // string attributes
  std::map<std::string, std::string> attStr;

  // int attributes
  std::map<std::string, int> attInt;

  // double attributes
  std::map<std::string, double> attDbl;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncattributes_new(NcAttributes_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncattributes_del(NcAttributes_t** self);


/**
 * Read the attributes of a variable
 * @param ncid netcdf id
 * @param varid variable id
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncattributes_read(NcAttributes_t** self, int ncid, int varid);


/**
 * Check if the variable is intensive
 * @return 1 if intensive, 0 otherwise
 */
extern "C"
int mnt_ncattributes_isIntensive(NcAttributes_t** self);

/**
 * Write the attributes of a variable
 * @param ncid netcdf id
 * @param varid variable id
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncattributes_write(NcAttributes_t** self, int ncid, int varid);

/**
 * Print the attributes
 * @return error code (0 is OK)
 */
extern "C"
int mnt_ncattributes_print(NcAttributes_t** self);


#endif // MNT_NC_ATTRIBUTES
