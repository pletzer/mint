#include "mntLIBRARY_API.h"
#include <vmtCellLocator.h>
#include <mntGrid.h>

#ifndef MNT_VMT_CELL_LOCATOR_CAPI
#define MNT_VMT_CELL_LOCATOR_CAPI

/**
 * C API wrapping vmtCellLocator (the abstract interface) and its two
 * concrete implementations, vmtLonLatCellLocator and vmtXYZCellLocator --
 * lets a caller (e.g. Python, via mint.LonLatCellLocator/mint.XYZCellLocator)
 * build and configure a locator of either kind directly, then hand it to
 * mnt_vectorinterp_setLocator (VectorInterp.setLocator) as a borrowed
 * reference -- VectorInterp does not take ownership or delete it, so the
 * caller must still call mnt_vmtcelllocator_del themselves once done.
 *
 * mnt_vectorinterp_buildLocator's useXYZLocator flag remains the simpler,
 * one-call way to get a locator VectorInterp both builds and owns; this
 * API is for when you want to build one yourself instead -- to configure it
 * beyond what buildLocator exposes, or to share one locator across several
 * VectorInterp/PolylineIntegral objects (see vmtCellLocator.h's docstring).
 *
 * All functions here operate through vmtCellLocator's shared (abstract)
 * interface, so the SAME set of calls (setDataSet, setNumberOfCellsPerBucket,
 * buildLocator, del) works regardless of which of the two `_new` variants
 * created the locator. setPeriodicityLengthX/enableFolding/setCubedSphere
 * are meaningful for a locator built via mnt_lonlatcelllocator_new only --
 * calling them on one built via mnt_xyzcelllocator_new logs an error rather
 * than silently doing the wrong thing (see vmtXYZCellLocator.h).
 */

/**
 * Constructor: a locator for a (lon, lat[, elev=0]) grid, periodic in
 * longitude and (optionally) folding across the poles.
 * @param self instance (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_lonlatcelllocator_new(vmtCellLocator** self);

/**
 * Constructor: a locator for a genuinely 3D-embedded (x, y, z) surface
 * mesh, no periodic seam -- see vmtXYZCellLocator.h.
 * @param self instance (output)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_xyzcelllocator_new(vmtCellLocator** self);

/**
 * Destructor -- works on a locator built by either constructor above.
 * @param self instance
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_del(vmtCellLocator** self);

/**
 * Attach the grid the locator will search.
 * @param self instance
 * @param grid grid (borrowed reference)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_setDataSet(vmtCellLocator** self, Grid_t* grid);

/**
 * Set average number of cells per bucket.
 * @param self instance
 * @param n number
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_setNumberOfCellsPerBucket(vmtCellLocator** self, int n);

/**
 * Set the periodicity length in x (only meaningful for a locator built via
 * mnt_lonlatcelllocator_new; logs an error, doesn't crash, for one built
 * via mnt_xyzcelllocator_new with periodX > 0 -- see vmtXYZCellLocator.h).
 * @param self instance
 * @param periodX length (0 means not periodic)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_setPeriodicityLengthX(vmtCellLocator** self, double periodX);

/**
 * Enable folding across poles (only meaningful for a locator built via
 * mnt_lonlatcelllocator_new).
 * @param self instance
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_enableFolding(vmtCellLocator** self);

/**
 * Declare whether the grid is a gnomonic cubed sphere (only meaningful for
 * a locator built via mnt_lonlatcelllocator_new).
 * @param self instance
 * @param isCubedSphere 1 if the grid is a cubed sphere, 0 otherwise
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_setCubedSphere(vmtCellLocator** self, int isCubedSphere);

/**
 * Build the locator -- call after setDataSet and any of the above.
 * @param self instance
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_vmtcelllocator_buildLocator(vmtCellLocator** self);

#endif // MNT_VMT_CELL_LOCATOR_CAPI
