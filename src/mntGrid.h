#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkGenericCell.h>
#include <vector>

#ifndef MNT_GRID
#define MNT_GRID

// index of each vertex
#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

struct Grid_t {

    const std::string EDGE_LENGTH_NAME = "_edge_lengths";

    // vertex raw data
    std::vector<double> verts;

    // VTK data needed to construct a grid
    vtkDoubleArray* pointData;
    vtkPoints* points;
    vtkUnstructuredGrid* grid;

    vtkUnstructuredGridReader* reader;

    // stores the fields (4 * numFaces)
    std::vector<vtkDoubleArray*> doubleArrays;

    // flat array of size numFaces * 4
    std::vector<size_t> faceNodeConnectivity;

    // flat array of size numFaces * 4
    std::vector<size_t> faceEdgeConnectivity;

    // flat array of size numEdges * 2
    std::vector<size_t> edgeNodeConnectivity;

    std::vector<double> edgeArcLengths;

    bool fixLonAcrossDateline;
    bool averageLonAtPole;

};


/**
 * Constructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C" 
int mnt_grid_new(Grid_t** self);

/**
 * Destructor
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_del(Grid_t** self);

/**
 * Set the grid flags
 * @param self instance of Grid_t
 * @param fixLonAcrossDateline set this to 1 if periodicty length should be added/subtracted to nodes in order to make the cell as compact as possible
 * @param averageLonAtPole set this to 1 if longitudes at the poles should take the average value of the node cell node's longitudes
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_setFlags(Grid_t** self, int fixLonAcrossDateline, int averageLonAtPole);

/**
 * Set the pointer to an array of points
 * @param self instance of Grid_t
 * @param nVertsPerCell number of vertices per cell
 * @param ncells number of cells
 * @param points flat array of size 4*ncells*3 (4 points per cell, each point has three coordinates)
 * @return error code (0 = OK)
 * @note the caller is responsible for managing the memory of the points array, which
 *       is expected to exist until the grid object is destroyed
 */
extern "C"
int mnt_grid_setPointsPtr(Grid_t** self, int nVertsPerCell, vtkIdType ncells, const double points[]);

/**
 * Attach data to the grid
 * @param self instance of Grid_t
 * @param varname data field name
 * @param nDataPerCell number of components per cell
 * @param points flat array of size ncells*nDataPerCell
 * @return error code (0 = OK)
 * @note the grid DOES NOT own the data - the caller is repsonsible for cleaning the data
 */
extern "C"
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]);

/**
 * Compute and store the arc length for each edge of the grid 
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_computeEdgeArcLengths(Grid_t** self);

/**
 * Get the VTK unstructured grid
 * @param self instance of Grid_t
 * @param grid_ptr address of vtlUnstructuredGrid object
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr);

/**
 * Load a grid from a 2D Ugrid file
 * @param self instance of Grid_t
 * @param fileAndMeshName column separated file and mesh name (e.g. cs_4.nc:physics)
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory when disposing of the grid
 * @note call mnt_grid_new prior to this call
 */
extern "C"
int mnt_grid_loadFrom2DUgrid(Grid_t** self, const char* fileAndMeshName);

/**
 * Load a grid from a VTK file
 * @param self instance of Grid_t
 * @param filename '\0' terminated file name
 * @return error code (0 = OK)
 * @note user should invoke mnt_grid_del to free memory when disposing of the grid
 * @note call mnt_grid_new prior to this call
 */
extern "C"
int mnt_grid_load(Grid_t** self, const char* filename);

/**
 * Dump the grid to a VTK file
 * @param self instance of Grid_t
 * @param filename '\0' terminated file name
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename);

/**
 * Print the content of the grid
 * @param self instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_print(Grid_t** self);

/**
 * Get the start or end point coordinates of an edge
 * @param self instance of Grid_t
 * @param cellId cell Id
 * @param edgeIndex edge index in the range 0-3
 * @param pointIndex point index in the range 0-1 (0 = start, 1 = end)
 * @param point the coordinates of the start point (output)
 * @param point the coordinates of the end point (output)
 * @return error code (0 = OK)
 */
extern "C"
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
extern "C" 
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
extern "C" 
int mnt_grid_getEdgeId(Grid_t** self, vtkIdType cellId, int edgeIndex, 
                       size_t* edgeId, int* signEdge);

/**
 * Get the nnumber of cells
 * @param self instance of Grid_t
 * @param numCells number of cells (output)
 * @return error code (0 = OK)
 */
extern "C"
int mnt_grid_getNumberOfCells(Grid_t** self, size_t* numCells);

/**
 * Get the number of edges
 * @param self instance of Grid_t
 * @param numEdges number of edges (output)
 * @return error code (0 = OK)
 * @note this will currently return 0 unless the grid was loaded from a UGRID 2D file
 */
extern "C"
int mnt_grid_getNumberOfEdges(Grid_t** self, size_t* numEdges);


#endif // MNT_GRID
