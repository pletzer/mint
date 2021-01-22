
#include <vector>
#include <map>
#include <string>
#include <vtkUnstructuredGrid.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>

#ifndef MNT_Polyline_INTEGRAL
#define MNT_Polyline_INTEGRAL

/**
 * @brief Edge-centred field regridding
 *
 * A class to compute flux integrals along a broken (poly-) line 
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
 * Constructs a polyline integral object
 * @param self instance of the polyline integral object
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
 * @param self instance of the polyline integral object'
 * @param grid instance of Grid_t
 * @return error code (0 = OK)
 */
extern "C"
int mnt_polylineintegral_setGrid(PolylineIntegral_t** self, Grid_t* grid);

/**
 * Set the target polyline 
 * @param self instance of the polyline integral object
 * @param npoints number of points
 * @param xyz flat array size npoints * 3 containing the xyz coordinates
 * @return error code (0 = OK)
 */
extern "C"
int mnt_polylineintegral_setPolyline(PolylineIntegral_t** self, int npoints, const double xyz[]);

/** 
 * Build the interpolation
 * @param self instance of the polyline integral object
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_build(PolylineIntegral_t** self);

/** 
 * Get the integral of the field along the line
 * @param self instance of the polyline integral object
 * @param data the edge integrated field values. 
 *        This array is expected to be dimensioned ncells * 4 (4 edges per cell). Each data value
 *        is a scalar representing the integral of the field over the edge. The directions of the 
 *        edge are (0, 0) -> (1, 0), (1, 0) -> (1, 1), (0, 1) -> (1, 1) and (0, 0) -> (0, 1) in
 *        parametric space
 * @param result line integral of the field over the polyline
 * @return error code (0 is OK)
 */
extern "C"
int mnt_polylineintegral_getIntegral(PolylineIntegral_t** self, const double data[], double* result);


#endif // MNT_POLY_LINE_INTEGRAL
