#include <vector>
#include <string>

#ifndef MNT_LATLON
#define MNT_LATLON

/**
 * A class to generate a uniform lat-lon grid conforming to the output of UM
 */

struct LatLon_t {
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
int mnt_latlon_new(LatLon_t** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_del(LatLon_t** self);

/**
 * Set the number of latitude cells
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLatCells(LatLon_t** self, size_t n);

/**
 * Set the number of longitude cells
 * @param n number of cells
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_setNumberOfLonCells(LatLon_t** self, size_t n);

/**
 * Build the grid
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_build(LatLon_t** self);

/**
 * Load a grid from file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_load(LatLon_t** self, const std::string& filename);

/**
 * Dump the grid to file
 * @param filename file name
 * @return error code (0 is OK)
 */
extern "C"
int mnt_latlon_dump(LatLon_t** self, const std::string& filename);

#endif // MNT_LATLON
