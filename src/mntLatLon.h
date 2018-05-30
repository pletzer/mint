#include <vector>
#include <string>

/**
 * A class to generate a uniform lat-lon grid conforming to the output of UM
 */

struct mntLatLon_t {
    std::vector<double> lats;
    std::vector<double> lons;
    double fillValue;
    double dLat;
    double dLon;
};

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_new(mntLatLon_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_del(mntLatLon_t** self);

/**
 * Set the number of latitude cells
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLatCells(mntLatLon_t** self, size_t n);

/**
 * Set the number of longitude cells
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLonCells(mntLatLon_t** self, size_t n);

/**
 * Build the grid
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_build(mntLatLon_t** self);

/**
 * Load a grid from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_load(mntLatLon_t** self, const std::string& filename);

/**
 * Dump the grid to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_dump(mntLatLon_t** self, const std::string& filename);


