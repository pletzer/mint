#include <mntPolylineIntegral.h>
#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>
#include <cmath>
#include <iostream>

double potential(const double p[]) {
    return p[1];
}

void testCartesian(size_t nx, size_t ny, double (*potentialFunc)(const double p[])) {

    int ier;

    Grid_t* grd = NULL;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);

    std::vector<double> verts(4*nx*ny*3); // number of cells times 4 nodes * 3 coordinates
    std::vector<double> data(4*nx*ny); // number of cells times 4 edges
    double dx = 1.0/double(nx);
    double dy = 1.0/double(ny);
    size_t k = 0;
    for (size_t i = 0; i < nx; ++i) {

        double x0 = 0. + dx*i;
        double x1 = x0 + dx;

        for (size_t j = 0; j < ny; ++j) {

            double y0 = 0. + dy*j;
            double y1 = y0 + dy;

            verts[4*k + 0 ] = x0;
            verts[4*k + 1 ] = y0;
            verts[4*k + 2 ] = 0.;

            verts[4*k + 3 ] = x1;
            verts[4*k + 4 ] = y0;
            verts[4*k + 5 ] = 0.;

            verts[4*k + 6 ] = x1;
            verts[4*k + 7 ] = y1;
            verts[4*k + 8 ] = 0.;

            verts[4*k + 9 ] = x0;
            verts[4*k + 10] = y1;
            verts[4*k + 11] = 0.;

            data[4*k + 0] = potentialFunc((const double*) &verts[4*k + 3 ]) - potentialFunc((const double*) &verts[4*k + 0 ]);
            data[4*k + 1] = potentialFunc((const double*) &verts[4*k + 6 ]) - potentialFunc((const double*) &verts[4*k + 3 ]);
            data[4*k + 2] = potentialFunc((const double*) &verts[4*k + 6 ]) - potentialFunc((const double*) &verts[4*k + 9 ]);
            data[4*k + 3] = potentialFunc((const double*) &verts[4*k + 9 ]) - potentialFunc((const double*) &verts[4*k + 0 ]);

            k++;
        }
    }
    vtkIdType ncells = nx * ny;
    ier = mnt_grid_setPointsPtr(&grd, 4, ncells, &verts[0]);
    assert(ier == 0);

    PolylineIntegral_t* pli = NULL;
    ier = mnt_polylineintegral_new(&pli);
    assert(ier == 0);

    ier = mnt_polylineintegral_setGrid(&pli, grd);
    assert(ier == 0);

    const int npoints = 4;
    std::vector<double> xyz({0., 0., 0.,
                             1., 0., 0.,
                             1., 1., 0.,
                             0., 1., 0.});
    ier = mnt_polylineintegral_setPolyline(&pli, npoints, (const double*) &xyz[0]);
    assert(ier == 0);

    ier = mnt_polylineintegral_build(&pli);
    assert(ier == 0);

    double totalFlux;
    ier = mnt_polylineintegral_getIntegral(&pli, (const double*) &data[0], &totalFlux);
    double totalFluxExact = potential(&xyz[0*3]) - potential(&xyz[0*3]);
    std::cout << "testCartesian: total flux = " << totalFlux << " exact: " << totalFluxExact << 
                 " error: " << totalFlux - totalFluxExact << '\n';
    assert(std::abs(totalFlux - totalFluxExact) < 1.e-10);

    ier = mnt_polylineintegral_del(&pli);
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0);
}

int main(int argc, char** argv) {

    testCartesian(10, 11, potential);

    return 0;
}
