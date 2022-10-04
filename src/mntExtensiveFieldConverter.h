#include "mntLIBRARY_API.h"
#include <mntGlobal.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVecN.h>
#include <sstream> // std::stringstream
#include "mntLogger.h"

#ifndef MNT_EXTENSIVE_FIELD_CONVERTER
#define MNT_EXTENSIVE_FIELD_CONVERTER


struct ExtensiveFieldConverter_t {
    // grid
    Grid_t* grid;
};

/**
 * Constructor
 * @param self instance of ExtensiveFieldConverter_t
 * @return error code (0 = OK)
 */
LIBRARY_API 
int mnt_extensivefieldconverter_new(ExtensiveFieldConverter_t** self);

/**
 * Destructor
 * @param self instance of ExtensiveFieldConverter_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldconverter_del(ExtensiveFieldConverter_t** self);

/**
 * Set the grid
 * @param self instance of ExtensiveFieldConverter_t
 * @param grid grid (borrowed reference)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldconverter_setGrid(ExtensiveFieldConverter_t** self, Grid_t* grid);

/**
 * Get the edge integrated data
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on edges, array of size numEdges
 * @param vy y component of vectors on edges, array of size numEdges
 * @param placement MNT_CELL_BY_CELL_DATA if data are cell by cell (size num cells * MNT_NUM_EDGES_PER_QUAD),
 *                  data are assumed to be on unique edges 
 *                  otherwise (size num edges)
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getEdgeData(ExtensiveFieldConverter_t** self,
                                            const double vx[], const double vy[], int placement,
                                            double data[]);

/**
 * Get the face integrated data
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on faces, array of size numEdges
 * @param vy y component of vectors on faces, array of size numEdges
 * @param placement MNT_CELL_BY_CELL_DATA if data are cell by cell (size num cells * MNT_NUM_EDGES_PER_QUAD),
 *                  data are assumed to be on unique edges 
 *                  otherwise (size num edges)
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getFaceData(ExtensiveFieldConverter_t** self,
                                            const double vx[], const double vy[], int placement,
                                            double data[]);

/**
 * Get the edge cell by cell data from the unique edge (or face) data
 * @param self instance of ExtensiveFieldConverter_t
 * @param uniqueEdgeData, input array of size numEdges
 * @param data cell by cell data, to be filled out
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldconverter_getEdgeCellByCellDataFromUniqueEdgeData(ExtensiveFieldConverter_t** self,
                                            const double uniqueEdgeData[],
                                            double data[]);

/**
 * Get the face cell by cell data from the unique edge (or face) data
 * @param self instance of ExtensiveFieldConverter_t
 * @param uniqueEdgeData, input array of size numEdges
 * @param data cell by cell data, to be filled out
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldconverter_getFaceCellByCellDataFromUniqueEdgeData(ExtensiveFieldConverter_t** self,
                                            const double uniqueEdgeData[],
                                            double data[]);


LIBRARY_API
int mnt_extensivefieldconverter__getEdgeDataFromCellByCellVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

LIBRARY_API
int mnt_extensivefieldconverter__getFaceDataFromCellByCellVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

LIBRARY_API
int mnt_extensivefieldconverter__getEdgeDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

LIBRARY_API
int mnt_extensivefieldconverter__getFaceDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

#endif // MNT_EXTENSIVE_FIELD_CONVERTER
