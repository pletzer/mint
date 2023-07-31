#include "mntLIBRARY_API.h"
#include <limits> // required by vtkUnstructuredGrid
#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>

#ifndef MNT_POLY_LINE_INTEGRAL
#define MNT_POLY_LINE_INTEGRAL

/**
 * A class to compute flux and line integrals along a broken (poly-) line. 
 */

struct PolylineIntegral_t {

    /** interpolation weights map: (cellId, edgeIndex) -> weight */
    std::map< std::pair<vtkIdType, int>, double > weights;

    /* grid */
    Grid_t* gridObj;

    /* cell locator */
    vmtCellLocator* loc;

    /* VTK grid (borrowed pointer) */
    vtkUnstructuredGrid* vgrid;
};

/**
 * Constructs a polyline integral object
 * @param self instance of the polyline integral object
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_new(PolylineIntegral_t** self);

/**
 * Destructor
 * @param self this instance
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_del(PolylineIntegral_t** self);

/** 
 * Set the grid
 * @param self instance of the polyline integral object
 * @param grid instance of Grid_t
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_setGrid(PolylineIntegral_t** self, Grid_t* grid);


/** 
 * Build the locator
 * @param self instance of the polyline integral object
 * @param numCellsPerBucket number of cells per bucket (small number allows a faster search,
 *                          large numbers are sometimes required, we recommend ~100)
 * @param periodX periodicity length (0. if not periodic)
 * @param enableFolding whether (1) or not (0) to enable folding across the poles
 *                      (allow |latitude| > 90)
 * @return error code (0 is OK)
 * @note call this after mnt_polylineintegral_setGrid
 */
LIBRARY_API
int mnt_polylineintegral_buildLocator(PolylineIntegral_t** self,
                                      int numCellsPerBucket,
                                      double periodX,
                                      int enableFolding);

/** 
 * Compute the interpolation weights
 * @param self instance of the polyline integral object
 * @param npoints number of points
 * @param xyz flat array size npoints * 3 containing the xyz coordinates
 * @param counterclock 1 if the edges go counterclockwise (that is [0,0]->[1,0], 
 *        [1,0]->[1,1], [1,1]->[0,1] and [0,1]->[0,0,]). Set counterclock = 0 if 
 *        the edges are oriented [0,0,]->[1,0], [1,0], [1,1], [0, 1]->[1,1] and 
 *        [0,0]->[0,1]
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_computeWeights(PolylineIntegral_t** self,
                                        int npoints, const double xyz[], int counterclock);

/** 
 * Get the integral of the field along the line
 * @param self instance of the polyline integral object
 * @param data edge integrated (extensive) field values.
 * @param placement MNT_CELL_BY_CELL_DATA if data are cell by cell 
 *                  (size of array is numCells * MNT_NUM_EDGES_PER_QUAD),
 *                  otherwise data are assumed to be on unique edges 
 *                  (size is numEdges)
 * @param result line integral of the field over the polyline
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self,
                                     const double data[], int placement,
                                     double* result);

/** 
 * Get the integral of the field along the line from a vector field
 * @param self instance of the polyline integral object
 * @param u eastward velocity on edges, size is num edges
 * @param v northward velocity on edges, size is num edges
 * @param fs function space (MNT_FUNC_SPACE_W1 or MNT_FUNC_SPACE_W2)
 * @param result line integral of the field over the polyline
 * @return error code (0 is OK)
 */
LIBRARY_API
int mnt_polylineintegral_vectorGetIntegral(PolylineIntegral_t** self,
                                           const double u[], const double v[], int fs, double* result);

#endif // MNT_POLY_LINE_INTEGRAL
