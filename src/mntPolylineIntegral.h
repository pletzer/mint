
#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntQuadEdgeIter.h>
#include <mntNcAttributes.h>
#include <mntNcFieldRead.h>
#include <mntNcFieldWrite.h>
#include <mntMultiArrayIter.h>

#ifndef MNT_Polyline_INTEGRAL
#define MNT_Polyline_INTEGRAL

/**
 * @brief Edge-centred field regridding
 *
 * A class to compute the regridding weights of an edge-centred field
 */

struct PolylineIntegral_t {

    /** pointer to the grid */
    Grid_t* grid;

    /** line coordinates */
    std::vector<double> lineXYZ;

    /** interpolation weights map: (cellId, edgeIndex) -> weight */
    std::map< std::pair<vtkIdType, int>, double > weights;

};

/**
 * Constructs a regridding object for edge centred fields
 * @param self instance of the regridding object
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_new(PolylineIntegral_t** self);

/**
 * Destructor
 * @param self this instance
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_del(PolylineIntegral_t** self);

/**
 * Set the grid
 * @param self instance of the regridding object'
 * @param grid instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_polylineintegral_setGrid(PolylineIntegral_t** self, Grid_t* grid);

/**
 * Set the target Polyline 
 * @param self instance of the regridding object
 * @param npoints number of points
 * @param xyz flat array of xyz coordinates
 * @return error code (0 = OK)
 */
extern "C"
int mnt_polylineintegral_setPolyline(PolylineIntegral_t** self, int npoints, const double xyz[]);

/** 
 * Build the interpolation
 * @param self instance of the regridding object
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_build(PolylineIntegral_t** self);

/** 
 * Get the integral of the field along the line
 * @param self instance of the regridding object
 * @param data the field values integrated over the grid's edges
 * @param result line integral of the filled over the Polyline
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self, const double data[], double* result);


#endif // MNT_POLY_LINE_INTEGRAL
