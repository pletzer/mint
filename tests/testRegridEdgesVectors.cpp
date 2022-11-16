#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntRegridEdges.h>
#include <mntLogger.h>
#include "saveEdgeVectors.h"
#include "saveEdgeVectorsXYZ.h"

#undef NDEBUG // turn on asserts
#include <cassert>

void test1() {

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
            double lam0 = p0[0] * deg2rad;
            double lam1 = p1[0] * deg2rad;
            double the0 = p0[1] * deg2rad;
            double the1 = p1[1] * deg2rad;

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

            // // rescale
            // dst_u[edgeId] *= rad2deg;
            // dst_v[edgeId] *= rad2deg;

            error += std::fabs(dst_u[edgeId] - dst_uExact) + std::fabs(dst_v[edgeId] - dst_vExact);
            // std::cerr << icell << ' ' << ie << " u=" << dst_u[edgeId] << " (" << dst_uExact << ") v=" << dst_v[edgeId] << " (" << dst_vExact << ")\n";
        }
    }
    error /= (dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    std::cerr << "test1: error = " << error << '\n';

    saveEdgeVectors(src_grd, src_u, src_v, "test1_src_vectors.vtk");
    saveEdgeVectors(dst_grd, dst_u, dst_v, "test1_dst_vectors.vtk");

    saveEdgeVectors(src_grd, src_u, src_v, "test1_src_vectorsXYZ.vtk");
    saveEdgeVectors(dst_grd, dst_u, dst_v, "test1_dst_vectorsXYZ.vtk");

    assert(error < 0.025);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}


void test2() {

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

            // stream function
            double s0 = cos(the0) * cos(lam0);
            double s1 = cos(the1) * cos(lam1);
            double dstDataExact = edgeSign * (s1 - s0);

            error += std::fabs(dstData[edgeId] - dstDataExact);
            // std::cerr << icell << ' ' << ie << " data=" << dstData[edgeId] << " (" << dstDataExact << ")\n";
        }
    }
    error /= (dst_numCells * MNT_NUM_EDGES_PER_QUAD);
    std::cerr << "test2: error = " << error << '\n';
    assert(error < 0.00015);

    // cleanup
    mnt_regridedges_del(&rgd);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}

int main(int argc, char** argv) {

    test2();
    test1();

    return 0;
}
