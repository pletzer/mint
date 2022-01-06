#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <mntGrid.h>
#include <mntNcFieldRead.h>
#include <mntVectorInterp.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <vector>


std::vector<double> computeFluxes(Grid_t* grd, 
                                  const std::vector<double>& u,
                                  const std::vector<double>& v) {


    std::size_t numCells, edgeId;
    int signEdge, ier;
    Vec3 p0, p1;

    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(ier == 0);

    std::vector<double> fluxes(numCells * 4);

    for (auto icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &signEdge);
            assert(ier == 0);

            ier = mnt_grid_getPoints(&grd, icell, ie, &p0[0], &p1[0]);
            assert(ier == 0);

            double dlon = p1[0] - p0[0];
            double dlat = p1[1] - p0[1];

            fluxes[4*icell + ie] = u[edgeId]*dlat - v[edgeId]*dlon;
        }
    }

    return fluxes;
}


void createGridAndLocator(Grid_t** grd, vmtCellLocator** cloc) {

    int ier;
    ier = mnt_grid_new(grd);
    assert(ier == 0);

    int fixLonAcrossDateline = 1;
    int averageLonAtPole = 1;
    int degrees = 1;
    ier = mnt_grid_setFlags(grd, fixLonAcrossDateline, averageLonAtPole, degrees);
    assert(ier == 0);

    // filename$meshname
    ier = mnt_grid_loadFromUgrid2D(grd, "${CMAKE_SOURCE_DIR}/data/lfric_diag_wind.nc$Mesh2d");
    assert(ier == 0);

    std::size_t numBadCells = 0;
    ier = mnt_grid_check(grd, &numBadCells);
    assert(numBadCells == 0);

    vtkUnstructuredGrid* ugrid;
    ier = mnt_grid_get(grd, &ugrid);
    assert(ier == 0);

    *cloc = vmtCellLocator::New();
    (*cloc)->SetDataSet(ugrid);
    (*cloc)->SetNumberOfCellsPerBucket(100);
    (*cloc)->BuildLocator();
}

void getVectors(Grid_t* grd, vmtCellLocator* cloc, const std::vector<double>& fluxes) {

    int ier;

    VectorInterp_t* vp = NULL;
    ier = mnt_vectorinterp_new(&vp);
    assert(ier == 0);
    ier = mnt_vectorinterp_setGrid(&vp, grd);
    assert(ier == 0);
    ier = mnt_vectorinterp_setLocator(&vp, cloc);
    assert(ier == 0);
    const double tol2 = 1.e-10;

    Vec3 p0, p1, pm, velocity;

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(ier == 0);

    // compute the velocity from fluxes at the edges and compare
    for (vtkIdType icell = 0; icell < numCells; ++icell) {

        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getPoints(&grd, icell, ie, &p0[0], &p1[0]);
            assert(ier == 0);

            // mid edge location
            pm = 0.5*(p0 + p1);

            // find the location of the point
            ier = mnt_vectorinterp_findPoints(&vp, 1, &pm[0], tol2);
            assert(ier == 0);

            // estimate the vector at the edge location
            ier = mnt_vectorinterp_getFaceVectorsFromCellByCellData(&vp, &fluxes[0],
                                                            &velocity[0]);
            assert(ier == 0);

            std::cout << "cell " << icell << " edge " << ie << " p0;p1=" << p0 << ";" << p1 << " flx=" << fluxes[4*icell+ie] << " v=" << velocity << '\n';

        }
    }

    ier = mnt_vectorinterp_del(&vp);
    assert(ier == 0);

}


void testZonal() {

    Grid_t* grd;
    vmtCellLocator* cloc;
    int ier;

    createGridAndLocator(&grd, &cloc);

    std::size_t numEdgeIds;
    ier = mnt_grid_getNumberOfEdges(&grd, &numEdgeIds);
    assert(ier == 0);

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(ier == 0);

    // zonal and meridional components of the velocity
    std::vector<double> u(numEdgeIds);
    std::vector<double> v(numEdgeIds);
    std::size_t edgeId;
    int signEdge;


    // set the velocity
    for (vtkIdType icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &signEdge);
            assert(ier == 0);

            u[edgeId] = 1.0; // 1 deg / time
            v[edgeId] = 0.0;
        }
    }

    // compute the fluxes
    std::vector<double> fluxes = computeFluxes(grd, u, v); // result is cell by cell 

    getVectors(grd, cloc, fluxes);


    // clean up
    cloc->Delete();

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}

void testMeridional() {

    Grid_t* grd;
    vmtCellLocator* cloc;
    int ier;

    createGridAndLocator(&grd, &cloc);

    std::size_t numEdgeIds;
    ier = mnt_grid_getNumberOfEdges(&grd, &numEdgeIds);
    assert(ier == 0);

    std::size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(ier == 0);

    // zonal and meridional components of the velocity
    std::vector<double> u(numEdgeIds);
    std::vector<double> v(numEdgeIds);
    std::size_t edgeId;
    int signEdge;


    // set the velocity
    for (auto icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < 4; ++ie) {

            ier = mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &signEdge);
            assert(ier == 0);

            u[edgeId] = 0.0; 
            v[edgeId] = 1.0; // 1 deg / time
        }
    }

    // compute the fluxes
    std::vector<double> fluxes = computeFluxes(grd, u, v); // result is cell by cell 

    getVectors(grd, cloc, fluxes);


    // clean up
    cloc->Delete();

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    testZonal();
    testMeridional();

    mnt_printLogMessages();
    return 0;
}
