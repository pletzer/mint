#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntExtensiveFieldAdaptor.h>
#include <mntVectorInterp.h>
#include <mntLogger.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void testFromUniqueEdgesCartesian() {

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
    std::vector<double> data(4);

    ier = mnt_extensivefieldadaptor_fromVectorField(&efa, &u[0], &v[0], &data[0], MNT_UNIQUE_EDGE_DATA, MNT_FUNC_SPACE_W1);
    assert(ier == 0);


    for (auto edgeId = 0; edgeId < 4; ++edgeId) {
        std::cout << "testFromUniqueEdgesCartesian: edge=" << edgeId << " extensive val=" << data[edgeId] << "\n";
    }

    const double tol = 1.e-15;

    // check
    assert(fabs(data[0] - (-10.0)) < tol);
    assert(fabs(data[1] - (-2.0)) < tol);
    assert(fabs(data[2] - (+30.0)) < tol);
    assert(fabs(data[3] - (-4.0)) < tol);

    mnt_grid_del(&grd);
    mnt_extensivefieldadaptor_del(&efa);
}



int main(int argc, char** argv) {

    testFromUniqueEdgesCartesian();

    mnt_printLogMessages();

    return 0;
}
