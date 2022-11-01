#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>
#include <mntGrid.h>
#include <mntExtensiveFieldAdaptor.h>
#include <mntRegridEdges.h>
#include <mntLogger.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void test2Cells() {

    /*

    x5  <2    x4   <6   x3
     
    3         4         ^
    v         v         5

    x0   <1   x1  0>    x2
    */

    std::vector<double> points({
        0., 0., 0.,
        2., 0., 0.,
        3., 0., 0.,
        3., 1., 0.,
        2., 1., 0.,
        0., 1., 0.
    });

    std::vector<std::size_t> face2nodes({
        0, 1, 4, 5,
        1, 2, 3, 4
    });
    std::vector<std::size_t> edge2nodes({
        1, 2,
        1, 0,
        4, 5,
        5, 0,
        4, 1,
        2, 3,
        3, 4
    });

    Grid_t* grd = NULL;
    mnt_grid_new(&grd);

    std::size_t numCells = 2;
    std::size_t numEdges = 7;
    std::size_t numPoints = 6;

    // Cartesian
    int ier = mnt_grid_setFlags(&grd, 0, 0, 0);
    assert(ier == 0);

    ier = mnt_grid_loadFromUgrid2DData(&grd, numCells, numEdges, numPoints, 
                                       &points[0], &face2nodes[0], &edge2nodes[0]);

     mnt_grid_print(&grd);

    ExtensiveFieldAdaptor_t* efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&efa);

    ier = mnt_extensivefieldadaptor_setGrid(&efa, grd);
    assert(ier == 0);

    std::vector<double> u({10., 11., 12., 13., 14., 15., 16.});
    std::vector<double> v({100., 110., 120., 130., 140., 150., 160.});
    std::vector<double> edgeData(numEdges); // W1
    std::vector<double> faceData(numEdges); // W2
    std::vector<double> u2(numEdges), v2(numEdges);

    // compute the extensive fields
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &edgeData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &faceData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    for (auto edgeId = 0; edgeId < numEdges; ++edgeId) {
        std::cout << "test2Cells: edge Id=" << edgeId << " edge =" << edgeData[edgeId] << 
                                                      " face =" << faceData[edgeId] << "\n";
    }

    const double tol = 1.e-15;

    // check edge integrals
    assert(fabs(edgeData[0] - (+10.)) < tol);
    assert(fabs(edgeData[1] - (-22.)) < tol);
    assert(fabs(edgeData[2] - (-24.)) < tol);
    assert(fabs(edgeData[3] - (-130.)) < tol);
    assert(fabs(edgeData[4] - (-140.)) < tol);
    assert(fabs(edgeData[5] - (+150.)) < tol);
    assert(fabs(edgeData[6] - (-16.)) < tol);

    // edge integrals
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &faceData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    for (auto edgeId = 0; edgeId < numEdges; ++edgeId) {
        std::cout << "test2Cells: edge Id=" << edgeId << " face extval=" << faceData[edgeId] << "\n";
    }

    // check face integrals
    assert(fabs(faceData[0] - (+100.)) < tol); // flux is positive but edge points down
    assert(fabs(faceData[1] - (-220.)) < tol);
    assert(fabs(faceData[2] - (-240.)) < tol);
    assert(fabs(faceData[3] - (-13.)) < tol);
    assert(fabs(faceData[4] - (-14.)) < tol);
    assert(fabs(faceData[5] - (+15.)) < tol);
    assert(fabs(faceData[6] - (-160.0)) < tol);

    // recover the vector field from the extensive field values
    ier = mnt_extensivefieldadaptor_toVectorField(&efa, &edgeData[0], &faceData[0],
                                                  &u2[0], &v2[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check
    for (auto i = 0; i < u.size(); ++i) {
        std::cout << "edge Id=" << i << " u = " << u[i] << " == " << u2[i] << 
                                     " v = " << v[i] << " == " << v2[i] << '\n';
        std::cout << " u error = " << fabs(u[i] - u2[i]) << " v error = " << fabs(v[i] - v2[i]) << '\n';
        assert(fabs(u[i] - u2[i]) < 1.e-7);
        assert(fabs(v[i] - v2[i]) < 1.e-7);
    }

    mnt_grid_del(&grd);
    mnt_extensivefieldadaptor_del(&efa);
}

void testCartesian() {

    /*

         3
    3  ..<.. 2
 0  V        ^ 2
    0  ..<.. 1
         1
    */

    std::vector<double> points({
        0., 0., 0.,
        1., 0., 0.,
        1., 1., 0.,
        0., 1., 0.,
    });

    std::vector<std::size_t> face2nodes({0, 1, 2, 3});
    std::vector<std::size_t> edge2nodes({3, 0,
                                         1, 0,
                                         1, 2,
                                         2, 3});

    Grid_t* grd = NULL;
    mnt_grid_new(&grd);

    std::size_t numCells = 1;
    std::size_t numEdges = 4;
    std::size_t numPoints = 4;

    // Cartesian
    int ier = mnt_grid_setFlags(&grd, 0, 0, 0);
    assert(ier == 0);

    ier = mnt_grid_loadFromUgrid2DData(&grd, numCells, numEdges, numPoints, 
                                       &points[0], &face2nodes[0], &edge2nodes[0]);

    ExtensiveFieldAdaptor_t* efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&efa);

    ier = mnt_extensivefieldadaptor_setGrid(&efa, grd);
    assert(ier == 0);


    std::vector<double> u({1., 2., 3., 4.});
    std::vector<double> v({10., 20., 30., 40.});
    std::vector<double> edgeData(numEdges); // W1
    std::vector<double> faceData(numEdges); // W2
    std::vector<double> u2(numEdges), v2(numEdges);

    // edge integrals
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &edgeData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);

    for (auto edgeId = 0; edgeId < numEdges; ++edgeId) {
        std::cout << "testCartesian: edge Id=" << edgeId << " edge extval=" << edgeData[edgeId] << "\n";
    }

    const double tol = 1.e-15;

    // check edge integrals
    assert(fabs(edgeData[0] - (-10.0)) < tol);
    assert(fabs(edgeData[1] - (-2.0)) < tol);
    assert(fabs(edgeData[2] - (+30.0)) < tol);
    assert(fabs(edgeData[3] - (-4.0)) < tol);

    // edge integrals
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &faceData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    for (auto edgeId = 0; edgeId < numEdges; ++edgeId) {
        std::cout << "testCartesian: edge Id=" << edgeId << " face extval=" << faceData[edgeId] << "\n";
    }

    // check face integrals
    assert(fabs(faceData[0] - (-1.0)) < tol); // flux is positive but edge points down
    assert(fabs(faceData[1] - (-20.0)) < tol); // flux is positive but edge points to the left
    assert(fabs(faceData[2] - (+3.0)) < tol);  // flux is positive and edge points up
    assert(fabs(faceData[3] - (-40.0)) < tol); // flux is positive but edge points to the left

    // recover the vector field from the extensive field values
    ier = mnt_extensivefieldadaptor_toVectorField(&efa, &edgeData[0], &faceData[0],
                                                  &u2[0], &v2[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check
    for (auto i = 0; i < u.size(); ++i) {
        std::cout << "edge=" << i << " u=" << u[i] << ' ' << u2[i] << '\n';
        std::cout << "edge=" << i << " v=" << v[i] << ' ' << v2[i] << '\n';
        assert(fabs(u[i] - u2[i]) < 1.e-8);
        assert(fabs(v[i] - v2[i]) < 1.e-8);
    }

    mnt_grid_del(&grd);
    mnt_extensivefieldadaptor_del(&efa);
}

void testLatLon2Itself() {

    int ier;

    // read the src/dst grids
    Grid_t* src_grd = NULL;
    mnt_grid_new(&src_grd);
    mnt_grid_setFlags(&src_grd, 0, 0, 1); // lat-lon, degrees
    mnt_grid_loadFromUgrid2DFile(&src_grd, "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc$latlon");
    std::size_t src_numEdges;
    mnt_grid_getNumberOfEdges(&src_grd, &src_numEdges);
    std::size_t src_numCells;
    mnt_grid_getNumberOfCells(&src_grd, &src_numCells);

    Grid_t* dst_grd = NULL;
    mnt_grid_new(&dst_grd);
    mnt_grid_setFlags(&dst_grd, 0, 0, 1); // lat-lon, degrees
    mnt_grid_loadFromUgrid2DFile(&dst_grd, "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc$latlon");
    std::size_t dst_numEdges;
    mnt_grid_getNumberOfEdges(&dst_grd, &dst_numEdges);
    std::size_t dst_numCells;
    mnt_grid_getNumberOfCells(&dst_grd, &dst_numCells);

    // set the vector field, which corresponds to streamfunction sin(theta) + cos(theta)*cos(lambda)
    std::vector<double> src_u(src_numEdges);
    std::vector<double> src_v(src_numEdges);
    std::size_t src_edgeId;
    Vec3 src_pt0, src_pt1;
    int src_edgeSign;
    for (vtkIdType src_cellId = 0; src_cellId < src_numCells; ++src_cellId) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&src_grd, src_cellId, ie, &src_edgeId, &src_edgeSign);
            mnt_grid_getPoints(&src_grd, src_cellId, ie, &src_pt0[0], &src_pt1[0]);
            Vec3 pm = 0.5*(src_pt0 + src_pt1);
            src_u[src_edgeId] = cos(pm[LAT_INDEX]*M_PI/180.);
            src_v[src_edgeId] = sin(pm[LON_INDEX]*M_PI/180.);
        }
    }

    std::vector<double> dst_u(dst_numEdges);
    std::vector<double> dst_v(dst_numEdges);
    std::size_t dst_edgeId;
    Vec3 dst_pt0, dst_pt1;
    int dst_edgeSign;
    for (vtkIdType dst_cellId = 0; dst_cellId < dst_numCells; ++dst_cellId) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&dst_grd, dst_cellId, ie, &dst_edgeId, &dst_edgeSign);
            mnt_grid_getPoints(&dst_grd, dst_cellId, ie, &dst_pt0[0], &dst_pt1[0]);
            Vec3 pm = 0.5*(dst_pt0 + dst_pt1);
            dst_u[dst_edgeId] = cos(pm[LAT_INDEX]*M_PI/180.);
            dst_v[dst_edgeId] = sin(pm[LON_INDEX]*M_PI/180.);
        }
    }

    std::vector<double> src_edgeIntegrals(src_numEdges);
    std::vector<double> src_faceIntegrals(src_numEdges);

    // compute edge/face integrals
    ExtensiveFieldAdaptor_t* src_efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&src_efa);

    ier = mnt_extensivefieldadaptor_setGrid(&src_efa, src_grd);
    assert(ier == 0);

    ier = mnt_extensivefieldadaptor_fromVectorField(&src_efa, &src_u[0], &src_v[0], &src_edgeIntegrals[0],
                                                    MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);
    ier = mnt_extensivefieldadaptor_fromVectorField(&src_efa, &src_u[0], &src_v[0], &src_faceIntegrals[0],
                                                    MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    // regrid
    std::vector<double> dst_edgeIntegrals(dst_numEdges);
    std::vector<double> dst_faceIntegrals(dst_numEdges);
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, src_grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dst_grd);
    assert(ier == 0);
    ier = mnt_regridedges_buildLocator(&rgd, 128, 360., 0);
    assert(ier == 0);
    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);
    ier = mnt_regridedges_apply(&rgd, &src_edgeIntegrals[0], &dst_edgeIntegrals[0],
        MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);
    ier = mnt_regridedges_apply(&rgd, &src_faceIntegrals[0], &dst_faceIntegrals[0],
        MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check that the src and dst extensive fields are the same
    double error = 0;
    for (auto i = 0; i < src_edgeIntegrals.size(); ++i) {
        std::cout << i << 
        " edge: " << src_edgeIntegrals[i] << ' ' << dst_edgeIntegrals[i] << 
        " face: " << src_faceIntegrals[i] << ' ' << dst_faceIntegrals[i] << '\n';
        error += fabs(src_edgeIntegrals[i] - dst_edgeIntegrals[i]) + 
                 fabs(src_faceIntegrals[i] - dst_faceIntegrals[i]);
    }
    assert(error < 1.e-8);

    // rebuild the vector field on the destination grid
    ExtensiveFieldAdaptor_t* dst_efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&dst_efa);

    ier = mnt_extensivefieldadaptor_setGrid(&dst_efa, dst_grd);
    assert(ier == 0);

    std::vector<double> dst_u2(dst_numEdges);
    std::vector<double> dst_v2(dst_numEdges);
    ier = mnt_extensivefieldadaptor_toVectorField(&dst_efa, 
                                                  &dst_edgeIntegrals[0], &dst_faceIntegrals[0],
                                                  &dst_u2[0], &dst_v2[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check that the rebuilt vector field is the same as the original one
    error = 0;
    for (auto i = 0; i < src_u.size(); ++i) {
        std::cout << i << 
        " u: " << src_u[i] << ' ' << dst_u2[i] << 
        " v: " << src_v[i] << ' ' << dst_v2[i] << '\n';
        error += fabs(src_u[i] - dst_u2[i]) + 
                 fabs(src_v[i] - dst_v2[i]);
    }
    assert(error < 1.e-8);


    mnt_extensivefieldadaptor_del(&dst_efa);
    mnt_regridedges_del(&rgd);
    mnt_extensivefieldadaptor_del(&src_efa);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}

void testCubedSphere2Itself() {

    int ier;

    // read the src/dst grids
    Grid_t* src_grd = NULL;
    mnt_grid_new(&src_grd);
    mnt_grid_setFlags(&src_grd, 1, 1, 1); // lat-lon, degrees
    mnt_grid_loadFromUgrid2DFile(&src_grd, "${CMAKE_SOURCE_DIR}/data/cs_4.nc$physics");
    std::size_t src_numEdges;
    mnt_grid_getNumberOfEdges(&src_grd, &src_numEdges);
    std::size_t src_numCells;
    mnt_grid_getNumberOfCells(&src_grd, &src_numCells);

    Grid_t* dst_grd = NULL;
    mnt_grid_new(&dst_grd);
    mnt_grid_setFlags(&dst_grd, 1, 1, 1); // cubed-sphere, degrees
    mnt_grid_loadFromUgrid2DFile(&dst_grd, "${CMAKE_SOURCE_DIR}/data/cs_4.nc$physics");
    std::size_t dst_numEdges;
    mnt_grid_getNumberOfEdges(&dst_grd, &dst_numEdges);
    std::size_t dst_numCells;
    mnt_grid_getNumberOfCells(&dst_grd, &dst_numCells);

    // set the vector field, which corresponds to streamfunction sin(theta) + cos(theta)*cos(lambda)
    std::vector<double> src_u(src_numEdges);
    std::vector<double> src_v(src_numEdges);
    std::size_t src_edgeId;
    Vec3 src_pt0, src_pt1;
    int src_edgeSign;
    for (vtkIdType src_cellId = 0; src_cellId < src_numCells; ++src_cellId) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&src_grd, src_cellId, ie, &src_edgeId, &src_edgeSign);
            mnt_grid_getPoints(&src_grd, src_cellId, ie, &src_pt0[0], &src_pt1[0]);
            Vec3 pm = 0.5*(src_pt0 + src_pt1);
            src_u[src_edgeId] = cos(pm[LAT_INDEX]*M_PI/180.);
            src_v[src_edgeId] = sin(pm[LON_INDEX]*M_PI/180.);
        }
    }

    std::vector<double> dst_u(dst_numEdges);
    std::vector<double> dst_v(dst_numEdges);
    std::size_t dst_edgeId;
    Vec3 dst_pt0, dst_pt1;
    int dst_edgeSign;
    for (vtkIdType dst_cellId = 0; dst_cellId < dst_numCells; ++dst_cellId) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&dst_grd, dst_cellId, ie, &dst_edgeId, &dst_edgeSign);
            mnt_grid_getPoints(&dst_grd, dst_cellId, ie, &dst_pt0[0], &dst_pt1[0]);
            Vec3 pm = 0.5*(dst_pt0 + dst_pt1);
            dst_u[dst_edgeId] = cos(pm[LAT_INDEX]*M_PI/180.);
            dst_v[dst_edgeId] = sin(pm[LON_INDEX]*M_PI/180.);
        }
    }

    std::vector<double> src_edgeIntegrals(src_numEdges);
    std::vector<double> src_faceIntegrals(src_numEdges);

    // compute edge/face integrals
    ExtensiveFieldAdaptor_t* src_efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&src_efa);

    ier = mnt_extensivefieldadaptor_setGrid(&src_efa, src_grd);
    assert(ier == 0);

    ier = mnt_extensivefieldadaptor_fromVectorField(&src_efa, &src_u[0], &src_v[0], &src_edgeIntegrals[0],
                                                    MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);
    ier = mnt_extensivefieldadaptor_fromVectorField(&src_efa, &src_u[0], &src_v[0], &src_faceIntegrals[0],
                                                    MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W2);
    assert(ier == 0);

    // regrid
    std::vector<double> dst_edgeIntegrals(dst_numEdges);
    std::vector<double> dst_faceIntegrals(dst_numEdges);
    RegridEdges_t* rgd = NULL;
    ier = mnt_regridedges_new(&rgd);
    assert(ier == 0);
    ier = mnt_regridedges_setSrcGrid(&rgd, src_grd);
    assert(ier == 0);
    ier = mnt_regridedges_setDstGrid(&rgd, dst_grd);
    assert(ier == 0);
    ier = mnt_regridedges_buildLocator(&rgd, 1024, 360., 0);
    assert(ier == 0);
    int debug = 2;
    ier = mnt_regridedges_computeWeights(&rgd, debug);
    assert(ier == 0);
    ier = mnt_regridedges_apply(&rgd, &src_edgeIntegrals[0], &dst_edgeIntegrals[0],
        MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);
    ier = mnt_regridedges_apply(&rgd, &src_faceIntegrals[0], &dst_faceIntegrals[0],
        MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check that the src and dst extensive fields are the same
    double error = 0;
    for (auto i = 0; i < src_edgeIntegrals.size(); ++i) {
        double e = fabs(src_edgeIntegrals[i] - dst_edgeIntegrals[i]) + 
                   fabs(src_faceIntegrals[i] - dst_faceIntegrals[i]);
        if (e > 1.e-8) {
            std::cout << "cs_4.nc$physics: " << i << 
            " edge: " << src_edgeIntegrals[i] << ' ' << dst_edgeIntegrals[i] << 
            " face: " << src_faceIntegrals[i] << ' ' << dst_faceIntegrals[i] << '\n';     
        }
        error += e;
    }
    std::cout << "error = " << error << '\n';
    assert(error < 1.e-6);

    // rebuild the vector field on the destination grid
    ExtensiveFieldAdaptor_t* dst_efa = NULL;
    ier = mnt_extensivefieldadaptor_new(&dst_efa);

    ier = mnt_extensivefieldadaptor_setGrid(&dst_efa, dst_grd);
    assert(ier == 0);

    std::vector<double> dst_u2(dst_numEdges);
    std::vector<double> dst_v2(dst_numEdges);
    ier = mnt_extensivefieldadaptor_toVectorField(&dst_efa, 
                                                  &dst_edgeIntegrals[0], &dst_faceIntegrals[0],
                                                  &dst_u2[0], &dst_v2[0], MNT_UNIQUE_EDGE_DATA);
    assert(ier == 0);

    // check that the rebuilt vector field is the same as the original one
    error = 0;
    for (auto i = 0; i < src_u.size(); ++i) {
        double e = fabs(src_u[i] - dst_u2[i]) + 
                 fabs(src_v[i] - dst_v2[i]);
        if (e > 1.e-8) {
            std::cout << "cs_4.nc$physics: " << i << 
            " u: " << src_u[i] << ' ' << dst_u2[i] << 
            " v: " << src_v[i] << ' ' << dst_v2[i] << '\n';
        }
        error += e;
    }
    std::cout << "error = " << error << '\n';
    assert(error < 1.e-8);


    mnt_extensivefieldadaptor_del(&dst_efa);
    mnt_regridedges_del(&rgd);
    mnt_extensivefieldadaptor_del(&src_efa);
    mnt_grid_del(&dst_grd);
    mnt_grid_del(&src_grd);
}


int main(int argc, char** argv) {

    testCubedSphere2Itself();
    // testLatLon2Itself();
    // testCartesian();
    // test2Cells();

    mnt_printLogMessages();

    return 0;
}
