#include "mntLIBRARY_API.h"
#include <limits> // required by vtkUnstructuredGrid
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkGenericCell.h>
#include <vector>
#include <mntGlobal.h>

#ifndef MNT_GRID
#define MNT_GRID

// index of each vertex
#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

struct Grid_t {

    const std::string EDGE_LENGTH_NAME = "_edge_lengths";

    // VTK data needed to construct a grid
    vtkDoubleArray* pointData;
    vtkPoints* points;
    vtkUnstructuredGrid* grid;

    vtkUnstructuredGridReader* reader;

    // stores the fields (MNT_NUM_EDGES_PER_QUAD * numFaces)
    std::vector<vtkDoubleArray*> doubleArrays;

    // flat array of size numFaces * MNT_NUM_VERTS_QUAD
    std::vector<std::size_t> faceNodeConnectivity;

    // flat array of size numFaces * MNT_NUM_EDGES_PER_QUAD
    std::vector<std::size_t> faceEdgeConnectivity;

    // flat array of size numEdges * MNT_NUM_VERTS_PER_EDGE
    std::vector<std::size_t> edgeNodeConnectivity;

    // edge arc lengths
    std::vector<double> edgeArcLengths;

    // number of unique edges, if loading from Ugrid2D
    std::size_t numEdges;

    // periodicity length, 360 if using degrees
    double periodX;

    // vertex raw data
    double *verts;

    // whether or not the vertices are owned by this instance
    bool ownsVerts;

    // whether we allow for a multiplicity of longitudes, useful 
    // for global domains
    bool fixLonAcrossDateline;

    // whether the longitudes should be adjusted at the poles to
    // minimize the cell sizes in lon-lat coordinates
    bool averageLonAtPole;
};

/**
 * Constructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_new(Grid_t** self);

/**
 * Destructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_del(Grid_t** self);

/**
 * Set the grid flags
 * @param self instance of Grid_t
 * @param fixLonAcrossDateline set this to 0 if periodicity length should NOT be added/subtracted to nodes in order to make the cell as compact as possible
 * @param averageLonAtPole set this to 0 if longitudes at the poles should NOT take the average value of the node cell node's longitudes
 * @param degrees set this to 0 if the coordinates are in radians
 * @note call this before reading the grid from file
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_setFlags(Grid_t** self, int fixLonAcrossDateline, int averageLonAtPole, int degrees);

/**
 * Set the pointer to an array of points
 * @param self instance of Grid_t
 * @param points flat array of size ncells*MNT_NUM_VERTS_PER_QUAD*3
 * @return error code (0 = OK)
 * @note the caller is responsible for managing the memory of the points array, which
 *       is expected to exist until the grid object is destroyed
 */
LIBRARY_API
int mnt_grid_setPointsPtr(Grid_t** self, double points[]);

/**
 * Get the pointer to the array of points
 * @param self instance of Grid_t
 * @param points pointer to the flat array of size ncells*MNT_NUM_VERTS_PER_QUAD*3
 *               holding the vertex coordinates
 * @return error code (0 = OK)
 * @see mnt_grid_setPointsPtr
 */
LIBRARY_API
int mnt_grid_getPointsPtr(Grid_t** self, double** points);


/**
 * Build the grid and connectivity
 * @param self instance of Grid_t
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_build(Grid_t** self, int nVertsPerCell, vtkIdType ncells);

/**
 * Attach data to the grid
 * @param self instance of Grid_t
 * @param varname data field name
 * @param nDataPerCell number of components per cell
 * @param data pointer to a flat array of size ncells*nDataPerCell
 * @return error code (0 = OK)
 * @note the grid DOES NOT own the data - the caller is responsible for ensuring that the data exist
 *       during the life of the grid object and for freeing the data
 */
LIBRARY_API
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]);

/**
 * Compute and store the arc length for each edge of the grid 
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 * @note assumes the sphere radius to be one
 */
LIBRARY_API
int mnt_grid_computeEdgeArcLengths(Grid_t** self);

/**
 * Get the length given a cell and an edge index
 * @param self instance of Grid_t
 * @param cellId cell Id
 * @param edgeIndex edge index in the range 0-3
 * @param res length of the edge assuming a radius of one
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_getEdgeArcLength(Grid_t** self, vtkIdType cellId, int edgeIndex, double* res);

/**
 * Get the VTK unstructured grid
 * @param self instance of Grid_t
 * @param grid_ptr address of vtlUnstructuredGrid object
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr);

/**
 * Load a grid from a 2D Ugrid file
 * @param self instance of Grid_t
 * @param fileAndMeshName column separated file and mesh name (e.g. cs_4.nc$physics)
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory when disposing of the grid
 * @note call mnt_grid_new prior to this call
 */
LIBRARY_API
int mnt_grid_loadFromUgrid2D(Grid_t** self, const char* fileAndMeshName);

/**
 * Load a grid from a VTK file
 * @param self instance of Grid_t
 * @param filename '\0' terminated file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory when disposing of the grid
 * @note call mnt_grid_new prior to this call
 */
LIBRARY_API
int mnt_grid_load(Grid_t** self, const char* filename);

/**
 * Dump the grid to a VTK file
 * @param self instance of Grid_t
 * @param filename '\0' terminated file name
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_dump(Grid_t** self, const char* filename);

/**
 * Print the content of the grid
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_print(Grid_t** self);

/**
 * Get the start or end point coordinates of an edge
 * @param self instance of Grid_t
 * @param cellId cell Id
 * @param edgeIndex edge index in the range 0-3
 * @param point0 the coordinates of the start point (output)
 * @param point1 the coordinates of the end point (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_getPoints(Grid_t** self, vtkIdType cellId, int edgeIndex, 
                       double point0[], double point1[]);

/**
 * Get the node Ids of a cellId, edge index pair
 * @param self instance of Grid_t
 * @param cellId cell Id
 * @param edgeIndex edge index in the range 0-3
 * @param nodeIds node Ids for the start/end vertices on the original grid (output)
 * @return error code (0 = OK)
 */
LIBRARY_API 
int mnt_grid_getNodeIds(Grid_t** self, vtkIdType cellId, int edgeIndex, vtkIdType nodeIds[]);

/**
 * Get the edge Id and direction of a cellId, edgeIndex pair
 * @param self instance of Grid_t
 * @param cellId cell Id
 * @param edgeIndex edge index in the range 0-3
 * @param edgeId edge Id (output)
 * @param signEdge +1 if edge edgeId points in the same direction as (cellid, edgeIndex), -1 otherwise
 * @return error code (0 = OK)
 */
LIBRARY_API 
int mnt_grid_getEdgeId(Grid_t** self, vtkIdType cellId, int edgeIndex, 
                       std::size_t* edgeId, int* signEdge);

/**
 * Get the number of cells
 * @param self instance of Grid_t
 * @param numCells number of cells (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_getNumberOfCells(Grid_t** self, std::size_t* numCells);

/**
 * Get the number of edges
 * @param self instance of Grid_t
 * @param numEdges number of edges (output)
 * @return error code (0 = OK)
 * @note this will currently return 0 unless the grid was loaded from a UGRID 2D file
 */
LIBRARY_API
int mnt_grid_getNumberOfEdges(Grid_t** self, std::size_t* numEdges);

/**
 * Check
 * @param self instance of Grid_t
 * @param numBadCells number of bad cells (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_grid_check(Grid_t** self, std::size_t* numBadCells);

#endif // MNT_GRID
