#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntExtensiveFieldAdaptor.h>
#include <mntVectorInterp.h>
#include <mntLogger.h>
#undef NDEBUG // turn on asserts
#include <cassert>

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
    std::vector<double> edgeData(4);
    std::vector<double> faceData(4);
    std::vector<double> u2(4), v2(4);

    // edge integrals
    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &edgeData[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);

    for (auto edgeId = 0; edgeId < 4; ++edgeId) {
        std::cout << "testFromUniqueEdgesCartesian: edge=" << edgeId << " edge extval=" << edgeData[edgeId] << "\n";
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

    for (auto edgeId = 0; edgeId < 4; ++edgeId) {
        std::cout << "testFromUniqueEdgesCartesian: edge=" << edgeId << " face extval=" << faceData[edgeId] << "\n";
    }

    // check face integrals
    assert(fabs(faceData[0] - (-1.0)) < tol); // flux is positive beut edge points down
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



int main(int argc, char** argv) {

    testCartesian();

    mnt_printLogMessages();

    return 0;
}
