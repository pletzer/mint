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
 * Get the edge integrated data from cell by cell edge vectors
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on edges, array of size numCells*MNT_NUM_EDGES_PER_QUAD
 * @param vy y component of vectors on edges, array of size numCells*MNT_NUM_EDGES_PER_QUAD
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getEdgeDataFromCellByCellVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

/**
 * Get the face integrated data from cell by cell face vectors
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on faces, array of size numCells*MNT_NUM_EDGES_PER_QUAD
 * @param vy y component of vectors on faces, array of size numCells*MNT_NUM_EDGES_PER_QUAD
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getFaceDataFromCellByCellVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

/**
 * Get the edge integrated data from unique edge vectors
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on edges, array of size numEdges
 * @param vy y component of vectors on edges, array of size numEdges
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getEdgeDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

/**
 * Get the face integrated data from unique edge vectors
 * @param self instance of ExtensiveFieldConverter_t
 * @param vx x component of vectors on faces, array of size numEdges
 * @param vy y component of vectors on faces, array of size numEdges
 * @param data edge integrated data, array size is numCells * MNT_NUM_EDGES_PER_QUAD
 * @return error code (0 = OK)
 * @note vx and vy should be compatible with the grid's coordinates. For instance, if vx and vy are 
 *       velocities in m/s and the coordinates are in degrees then one needs transform 
 *       vx -> vx*(180/pi)/(A*cos(lat*pi/180)) and vy -> vy*(180/pi)/A to get velocities in deg/sec
 *       (A is the planet radius).
 */
LIBRARY_API
int mnt_extensivefieldconverter_getFaceDataFromUniqueEdgeVectors(ExtensiveFieldConverter_t** self,
                                                              const double vx[], const double vy[],
                                                              double data[]);

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


#endif // MNT_EXTENSIVE_FIELD_CONVERTER