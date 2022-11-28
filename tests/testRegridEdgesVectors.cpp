#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntRegridEdges.h>
#include <mntLogger.h>
#include "saveEdgeVectors.h"
#include "saveEdgeVectorsXYZ.h"

#undef NDEBUG // turn on asserts
#include <cassert>

void testVector2Vector() {

    // testting vector to vector regridding

    int ier;

    // read the src/dst grids
    Grid_t* src_grd = NULL;
    ier = mnt_grid_new(&src_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&src_grd, 0, 0, 1); // lon-lat, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&src_grd, "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc$latlon");
    assert(ier == 0);

    Grid_t* dst_grd = NULL;
    ier = mnt_grid_new(&dst_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&dst_grd, 1, 1, 1); // cubed-sphere, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&dst_grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc$physics");
    assert(ier == 0);

    ier = mnt_grid_dump(&src_grd, "latlon100x50_grid.vtk");
    assert(ier == 0);
    ier = mnt_grid_dump(&dst_grd, "cs16_grid.vtk");
    assert(ier == 0);

    // set the u,v field on the src grid
    std::size_t src_numEdges;
    ier = mnt_grid_getNumberOfEdges(&src_grd, &src_numEdges);
    std::size_t src_numCells;
    ier = mnt_grid_getNumberOfCells(&src_grd, &src_numCells);
    std::cout << "src num cells: " << src_numCells << " src num edges: " << src_numEdges << '\n';

    std::size_t edgeId;
    int edgeSign;
    Vec3 p0, p1;
    std::vector<double> src_u(src_numEdges), src_v(src_numEdges);
    const double deg2rad = M_PI/180.;
    for (auto icell = 0; icell < src_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            ier = mnt_grid_getEdgeId(&src_grd, icell, ie, &edgeId, &edgeSign);
            ier = mnt_grid_getPoints(&src_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[LON_INDEX] * deg2rad;
            double lam1 = p1[LON_INDEX] * deg2rad;
            double the0 = p0[LAT_INDEX] * deg2rad;
            double the1 = p1[LAT_INDEX] * deg2rad;

            // // stream function
            // double s0 = cos(the0) * cos(lam0);
            // double s1 = cos(the1) * cos(lam1);

            // vector field
            double lamMid = 0.5*(lam0 + lam1);
            double theMid = 0.5*(the0 + the1);
            src_u[edgeId] = - sin(theMid) * cos(lamMid);
            src_v[edgeId] = sin(lamMid);
        }
    }

    // regridder
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, src_grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dst_grd);
    assert(ier == 0);
    int numCellsPerBucket = 256;
    double periodX = 360.;
    int enableFolding = 1;
    ier = mnt_regridedges_buildLocator(&rgd, numCellsPerBucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);

    std::size_t dst_numEdges;
    ier = mnt_grid_getNumberOfEdges(&dst_grd, &dst_numEdges);
    std::size_t dst_numCells;
    ier = mnt_grid_getNumberOfCells(&dst_grd, &dst_numCells);
    std::cout << "dst num cells: " << dst_numCells << " dst num edges: " << dst_numEdges << '\n';

    std::vector<double> dst_u(dst_numEdges), dst_v(dst_numEdges);
    ier = mnt_regridedges_vectorApply(&rgd, &src_u[0], &src_v[0],
                                            &dst_u[0], &dst_v[0],
                                            MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    // check
    const double rad2deg = 180./M_PI;
    double error = 0;
    for (auto icell = 0; icell < dst_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            ier = mnt_grid_getEdgeId(&dst_grd, icell, ie, &edgeId, &edgeSign);
            ier = mnt_grid_getPoints(&dst_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);

            // vector field
            double lamMid = 0.5*(lam0 + lam1);
            double theMid = 0.5*(the0 + the1);
            double dst_uExact = - sin(theMid) * cos(lamMid);
            double dst_vExact = sin(lamMid);

            error += std::fabs(dst_u[edgeId] - dst_uExact) + std::fabs(dst_v[edgeId] - dst_vExact);
            // std::cerr << icell << ' ' << ie << " u=" << dst_u[edgeId] << " (" << dst_uExact << ") v=" << dst_v[edgeId] << " (" << dst_vExact << ")\n";
        }
    }
    error /= (dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    std::cerr << "testVector2Vector: error = " << error << '\n';

    saveEdgeVectors(src_grd, src_u, src_v, "testVector2Vector_src_vectors.vtk");
    saveEdgeVectors(dst_grd, dst_u, dst_v, "testVector2Vector_dst_vectors.vtk");

    saveEdgeVectorsXYZ(src_grd, src_u, src_v, "testVector2Vector_src_vectorsXYZ.vtk");
    saveEdgeVectorsXYZ(dst_grd, dst_u, dst_v, "testVector2Vector_dst_vectorsXYZ.vtk");

    assert(error < 0.025);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}


void testExtensiveFieldRegriddingUniqueEdgeData() {

    // testing extensive field regridding using unique edge data

    int ier;

    // read the src/dst grids
    Grid_t* src_grd = NULL;
    ier = mnt_grid_new(&src_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&src_grd, 0, 0, 1); // lon-lat, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&src_grd, "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc$latlon");
    assert(ier == 0);

    Grid_t* dst_grd = NULL;
    ier = mnt_grid_new(&dst_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&dst_grd, 1, 1, 1); // cubed-sphere, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&dst_grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc$physics");
    assert(ier == 0);

    ier = mnt_grid_dump(&src_grd, "latlon100x50_grid.vtk");
    assert(ier == 0);
    ier = mnt_grid_dump(&dst_grd, "cs16_grid.vtk");
    assert(ier == 0);

    // set the data on the src grid
    std::size_t src_numEdges;
    ier = mnt_grid_getNumberOfEdges(&src_grd, &src_numEdges);
    std::size_t src_numCells;
    ier = mnt_grid_getNumberOfCells(&src_grd, &src_numCells);
    std::cout << "src num cells: " << src_numCells << " src num edges: " << src_numEdges << '\n';

    std::size_t edgeId;
    int edgeSign;
    Vec3 p0, p1;
    std::vector<double> srcData(src_numEdges);
    const double deg2rad = M_PI/180.;
    for (auto icell = 0; icell < src_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            ier = mnt_grid_getEdgeId(&src_grd, icell, ie, &edgeId, &edgeSign);
            ier = mnt_grid_getPoints(&src_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            srcData[edgeId] = edgeSign * (s1 - s0);

        }
    }

    // regridder
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, src_grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dst_grd);
    assert(ier == 0);
    int numCellsPerBucket = 256;
    double periodX = 360.;
    int enableFolding = 0;
    ier = mnt_regridedges_buildLocator(&rgd, numCellsPerBucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);

    std::size_t dst_numEdges;
    ier = mnt_grid_getNumberOfEdges(&dst_grd, &dst_numEdges);
    std::size_t dst_numCells;
    ier = mnt_grid_getNumberOfCells(&dst_grd, &dst_numCells);
    std::cout << "dst num cells: " << dst_numCells << " dst num edges: " << dst_numEdges << '\n';

    std::vector<double> dstData(dst_numEdges);
    std::vector<double> dstUFluxApprox(dst_numEdges);
    std::vector<double> dstVFluxApprox(dst_numEdges);
    std::vector<double> dstUFluxExact(dst_numEdges);
    std::vector<double> dstVFluxExact(dst_numEdges);
    ier = mnt_regridedges_apply(&rgd, &srcData[0], &dstData[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check
    double error = 0;
    for (auto icell = 0; icell < dst_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            ier = mnt_grid_getEdgeId(&dst_grd, icell, ie, &edgeId, &edgeSign);
            ier = mnt_grid_getPoints(&dst_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

            double lamMid = 0.5*(lam0 + lam1);
            double theMid = 0.5*(the0 + the1);

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            double dstDataExact = edgeSign * (s1 - s0);

            double uMid = - sin(theMid) * cos(lamMid);
            double vMid = sin(lamMid);
            dstUFluxApprox[edgeId] = uMid * (the1 - the0);
            dstVFluxApprox[edgeId] = - vMid * (lam1 - lam0) * cos(theMid);

            // int u d theta = - int sin(theta)*cos(lambda) dtheta
            dstUFluxExact[edgeId] = ((the0-the1)*((the0-the1)*cos(lam0)*cos(the0)-(the0-the1)*cos(lam1)*cos(the1)+(lam0-lam1)*(sin(lam0)*sin(the0)-sin(lam1)*sin(the1))))/((lam0-lam1+the0-the1)*(lam0-lam1-the0+the1));
            // - int  v cos(theta)*dlambda = - int sin(lambda)*cos(theta)* dlambda
            dstVFluxExact[edgeId] = -(((lam0-lam1)*((lam0-lam1)*cos(lam0)*cos(the0)+(-lam0+lam1)*cos(lam1)*cos(the1)+(the0-the1)*(sin(lam0)*sin(the0)-sin(lam1)*sin(the1))))/((lam0-lam1+the0-the1)*(lam0-lam1-the0+the1)));

            error += std::fabs(dstData[edgeId] - dstDataExact);
            // std::cerr << icell << ' ' << ie << " data=" << dstData[edgeId] << " (" << dstDataExact << ")\n";

            // DEBUG
            if (icell == 1160) {
                std::cout << "testExtensiveFieldRegriddingUniqueEdgeData: cell " << icell << " edge " << ie << " edge sign=" << edgeSign << " ext field=" << dstData[edgeId] << " (exact: " << dstDataExact << ") edge p0=" << p0 << " -> " << p1 << '\n';
            }

        }
    }
    error /= (dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    std::cerr << "testExtensiveFieldRegriddingUniqueEdgeData: error = " << error << '\n';
    saveEdgeVectors(dst_grd, dstUFluxApprox, dstVFluxApprox, "testExtensiveFieldRegriddingUniqueEdgeData_dst_approx_fluxes.vtk");
    saveEdgeVectors(dst_grd, dstUFluxExact, dstVFluxExact, "testExtensiveFieldRegriddingUniqueEdgeData_dst_exact_fluxes.vtk");

    assert(error < 0.00015);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}


void testExtensiveFieldCellByCellData() {

    // testing extensive field regridding using cell by cell data

    int ier;

    // read the src/dst grids
    Grid_t* src_grd = NULL;
    ier = mnt_grid_new(&src_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&src_grd, 0, 0, 1); // lon-lat, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&src_grd, "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc$latlon");
    assert(ier == 0);

    Grid_t* dst_grd = NULL;
    ier = mnt_grid_new(&dst_grd);
    assert(ier == 0);
    ier = mnt_grid_setFlags(&dst_grd, 1, 1, 1); // cubed-sphere, degrees
    assert(ier == 0);
    ier = mnt_grid_loadFromUgrid2DFile(&dst_grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc$physics");
    assert(ier == 0);


    std::size_t src_numCells;
    ier = mnt_grid_getNumberOfCells(&src_grd, &src_numCells);
    std::cout << "src num cells: " << src_numCells << '\n';

    Vec3 p0, p1;
    std::vector<double> srcData(src_numCells * MNT_NUM_EDGES_PER_QUAD);
    const double deg2rad = M_PI/180.;
    for (auto icell = 0; icell < src_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            ier = mnt_grid_getPoints(&src_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + ie;
            int sign = 1;
            srcData[k] = sign * (s1 - s0);

        }
    }

    // regridder
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, src_grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dst_grd);
    assert(ier == 0);
    int numCellsPerBucket = 256;
    double periodX = 360.;
    int enableFolding = 0;
    ier = mnt_regridedges_buildLocator(&rgd, numCellsPerBucket, periodX, enableFolding);
    assert(ier == 0);

    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);

    std::size_t dst_numCells;
    ier = mnt_grid_getNumberOfCells(&dst_grd, &dst_numCells);
    std::cout << "dst num cells: " << dst_numCells << '\n';

    std::vector<double> dstData(dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    ier = mnt_regridedges_apply(&rgd, &srcData[0], &dstData[0], MNT_CELL_BY_CELL_DATA);
    assert(ier == 0);

    // check
    double error = 0;
    for (auto icell = 0; icell < dst_numCells; ++icell) {
        for (auto ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

            ier = mnt_grid_getPoints(&dst_grd, icell, ie, &p0[0], &p1[0]);
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

            double lamMid = 0.5*(lam0 + lam1);
            double theMid = 0.5*(the0 + the1);

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            int sign = 1;
            if (ie == 0 || ie == 2) sign = -1;
            double dstDataExact = sign * (s1 - s0);

            std::size_t k = icell*MNT_NUM_EDGES_PER_QUAD + ie;
            error += std::fabs(dstData[k] - dstDataExact);
            // std::cerr << icell << ' ' << ie << " data=" << dstData[edgeId] << " (" << dstDataExact << ")\n";

            // DEBUG
            if (icell == 1160) {
                std::cout << "testExtensiveFieldCellByCellData: cell " << icell << " edge " << ie << " ext field=" << dstData[k] << " edge p0=" << p0 << " -> " << p1 << '\n';
            }

        }
    }
    error /= (dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    std::cerr << "testExtensiveFieldCellByCellData: error = " << error << '\n';
    assert(error < 0.05);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}

int main(int argc, char** argv) {

    testExtensiveFieldCellByCellData();
    testExtensiveFieldRegriddingUniqueEdgeData();
    testVector2Vector();

    return 0;
}
