#include "mntLIBRARY_API.h"
#include <mntGlobal.h>
#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVecN.h>
#include <sstream> // std::stringstream
#include "mntLogger.h"

#ifndef MNT_EXTENSIVE_FIELD_ADAPTOR
#define MNT_EXTENSIVE_FIELD_ADAPTOR


struct ExtensiveFieldAdaptor_t {

    // grid
    Grid_t* grid;

    // geometry
    int geo;

    std::size_t numCells;
    std::size_t numEdges;
};

/**
 * Constructor
 * @param self instance of extensivefieldadaptor_t
 * @return error code (0 = OK)
 */
LIBRARY_API 
int mnt_extensivefieldadaptor_new(ExtensiveFieldAdaptor_t** self);

/**
 * Destructor
 * @param self instance of extensivefieldadaptor_t
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldadaptor_del(ExtensiveFieldAdaptor_t** self);

/**
 * Set the grid
 * @param self instance of ExtensiveFieldAdaptor_t
 * @param grid grid (borrowed reference)
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldadaptor_setGrid(ExtensiveFieldAdaptor_t** self, Grid_t* grid);

/**
 * Get the extensive field from a vector field
 * @param self instance of ExtensiveFieldAdaptor_t
 * @param u x-component, size of array depends on placement, see below
 * @param v y-component, size of array depends on placement, see below
 * @param data extensive field, size of array depends on placement, see below (output)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. In the former case the size of the 
 *                  arrays should be num_cells*MNT_NUM_EDGES_PER_QUAD while in the latter it should be num_edges
 * @param fs function space, either MNT_FUNC_SPACE_W1 or MNT_FUNC_SPACE_W2
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldadaptor_fromVectorField(ExtensiveFieldAdaptor_t** self, const double u[], const double v[],
                                              double data[], int placement, int fs);

/**
 * Get the vector field from an extensive field array
 * @param self instance of ExtensiveFieldAdaptor_t
 * @param edgeData extensive edge field, size of array depends on placement, see below
 * @param faceData extensive face field, size of array depends on placement, see below
 * @param u x-component, size of array depends on placement, see below (output)
 * @param v y-component, size of array depends on placement, see below (output)
 * @param placement either MNT_CELL_BY_CELL_DATA or MNT_UNIQUE_EDGE_DATA. In the former case the size of the
 *                  arrays should be num_cells*MNT_NUM_EDGES_PER_QUAD while in the latter it should be num_edges
 * @return error code (0 = OK)
 */
LIBRARY_API
int mnt_extensivefieldadaptor_toVectorField(ExtensiveFieldAdaptor_t** self,
                                            const double edgeData[], const double faceData[],
                                            double u[], double v[],
                                            int placement);

// private methods
int mnt_extensivefieldadaptor__fromVectorFieldEdgeUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]);
int mnt_extensivefieldadaptor__fromVectorFieldFaceUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]);
int mnt_extensivefieldadaptor__fromVectorFieldEdgeCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]);
int mnt_extensivefieldadaptor__fromVectorFieldFaceCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                            const double u[], const double v[],
                                                            double data[]);
int mnt_extensivefieldadaptor__toVectorFieldUniqueIdData(ExtensiveFieldAdaptor_t** self,
                                                           const double edgeData[],
                                                           const double faceData[],
                                                           double u[], double v[]);
int mnt_extensivefieldadaptor__toVectorFieldCellByCellData(ExtensiveFieldAdaptor_t** self,
                                                           const double edgeData[],
                                                           const double faceData[],
                                                           double u[], double v[]);

#endif // MNT_EXTENSIVE_FIELD_ADAPTOR
