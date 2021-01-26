#include <mntPolylineIntegral.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

double potential(const double p[]) {
    const double x = p[0];
    const double y = p[1];
    return y*(2*x + 1);
}

void testCartesian(size_t nx, size_t ny, double (*potentialFunc)(const double p[])) {

    int ier;

    Grid_t* grd = NULL;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    std::vector<double> verts(4*nx*ny*3); // number of cells * 4 nodes * 3 coordinates
    std::vector<double> data(4*nx*ny); // number of cells * 4 edges
    double dx = 1.0/double(nx);
    double dy = 1.0/double(ny);

    // build the grid, which is made of independent quad cells. For each of the cells we
    // set the edge field values
    size_t k = 0;
    const double* p0;
    const double* p1;
    for (size_t i = 0; i < nx; ++i) {

        double x0 = 0. + dx*i;
        double x1 = x0 + dx;

        for (size_t j = 0; j < ny; ++j) {

            double y0 = 0. + dy*j;
            double y1 = y0 + dy;

            // set the vertices. Note 3d even when the quads are in the plane.
            verts[4*3*k + 0 ] = x0;
            verts[4*3*k + 1 ] = y0;
            verts[4*3*k + 2 ] = 0.;

            verts[4*3*k + 3 ] = x1;
            verts[4*3*k + 4 ] = y0;
            verts[4*3*k + 5 ] = 0.;

            verts[4*3*k + 6 ] = x1;
            verts[4*3*k + 7 ] = y1;
            verts[4*3*k + 8 ] = 0.;

            verts[4*3*k + 9 ] = x0;
            verts[4*3*k + 10] = y1;
            verts[4*3*k + 11] = 0.;

            // set the edge integrals, just a difference of potential
            // notre the orientation of the edges
            //      2
            //  3 -->---2
            //  |       |
            //  ^ 3     ^ 1
            //  |       |
            //  0 -->---1
            //      0

            // south edge
            p0 = (const double*) &verts[4*3*k + 3*0 ];
            p1 = (const double*) &verts[4*3*k + 3*1 ];
            data[4*k + 0] = potentialFunc(p1) - potentialFunc(p0);


            // east edge
            p0 = (const double*) &verts[4*3*k + 3*1 ];
            p1 = (const double*) &verts[4*3*k + 3*2 ];
            data[4*k + 1] = potentialFunc(p1) - potentialFunc(p0);

            // north edge
            p0 = (const double*) &verts[4*3*k + 3*3 ];
            p1 = (const double*) &verts[4*3*k + 3*2 ];
            data[4*k + 2] = potentialFunc(p1) - potentialFunc(p0);

            // west edge
            p0 = (const double*) &verts[4*3*k + 3*0 ];
            p1 = (const double*) &verts[4*3*k + 3*3 ];
            data[4*k + 3] = potentialFunc(p1) - potentialFunc(p0);

            // increment the cell index
            k++;
        }
    }
    vtkIdType ncells = nx * ny;
    ier = mnt_grid_setPointsPtr(&grd, 4, ncells, &verts[0]);
    assert(ier == 0);

    PolylineIntegral_t* pli = NULL;
    ier = mnt_polylineintegral_new(&pli);
    assert(ier == 0);

    // set the (open) contour to integrate the flux over
    const int npoints = 4;
    std::vector<double> xyz({0., 0., 0.,
                             1., 0., 0.,
                             1., 1., 0.,
                             0., 1., 0.});
    int counterclock = 0;
    ier = mnt_polylineintegral_build(&pli, grd, npoints, (const double*) &xyz[0], counterclock);
    assert(ier == 0);

    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0], &totalFlux);

    // exact flux is the difference of potential between end and starting points
    double totalFluxExact = potential(&xyz[(npoints - 1)*3]) - potential(&xyz[0*3]);

    std::cout << "testCartesian: total flux = " << totalFlux << " exact: " << totalFluxExact << 
                 " error: " << totalFlux - totalFluxExact << '\n';
    assert(std::abs(totalFlux - totalFluxExact) < 1.e-10);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    testCartesian(3, 2, potential);

    return 0;
}
