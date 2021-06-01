#include <vmtCellLocator.h>
#include <mntGrid.h>
#include <mntVectorInterp.h>
#undef NDEBUG // turn on asserts
#include <cassert>

/**
 * Create uniform VTK grid objects
 * @param nx number of x cells
 * @param ny number of y cells
 * @param points array of size (4*ny*nx) to be filled in
 */
void createUniformGrid(int nx, int ny, double points[]) {

    double dx = 360. / (double) nx;
    double dy = 180. / (double) ny;

    int k = 0;
    for (int i = 0; i < nx; ++i) {
        double x0 = i*dx;
        double x1 = x0 + dx;
        for (int j = 0; j < ny; ++j) {
            double y0 = -90.0 + j*dy;
            double y1 = y0 + dy;

            points[4*k*3 +  0] = x0;
            points[4*k*3 +  1] = y0;
            points[4*k*3 +  2] = 0.0;

            points[4*k*3 +  3] = x1;
            points[4*k*3 +  4] = y0;
            points[4*k*3 +  5] = 0.0;

            points[4*k*3 +  6] = x1;
            points[4*k*3 +  7] = y1;
            points[4*k*3 +  8] = 0.0;

            points[4*k*3 +  9] = x0;
            points[4*k*3 + 10] = y1;
            points[4*k*3 + 11] = 0.0;

            k++;
        }
    }
}

void testContainsPoint(int nx, int ny) {

    std::vector<double> points(4*ny*nx*3);
    createUniformGrid(nx, ny, &points[0]);

    Grid_t* mgrid = NULL;
    mnt_grid_new(&mgrid);
    int fixLonAcrossDateline = 0;
    int averageLonAtPole = 0;
    int degrees = 1;
    mnt_grid_setFlags(&mgrid, fixLonAcrossDateline, averageLonAtPole, degrees);
    mnt_grid_setPointsPtr(&mgrid, 4, ny*nx, &points[0]);

    // create locator
    vmtCellLocator* cloc = vmtCellLocator::New();
    vtkUnstructuredGrid* ugrid;
    mnt_grid_get(&mgrid, &ugrid);
    cloc->SetDataSet(ugrid);
    cloc->SetNumberOfCellsPerBucket(10);
    cloc->BuildLocator();
    cloc->setPeriodicityLengthX(360.0);

    // set some data
    std::vector<double> data(4*ny*nx);
    for (std::size_t i = 0; i < 4*ny*nx; ++i) {
        data[i] = (double) i;
    }

    Vec3 target, resVec;
    // inside the first face
    target[0] = 0.00001;
    target[1] = -89.99999;
    target[2] = 0.;

    VectorInterp_t* vp = NULL;
    mnt_vectorinterp_new(&vp);
    mnt_vectorinterp_set(&vp, mgrid, cloc);
    double tol2 = 1.e-10;
    mnt_vectorinterp_find(&vp, &target[0], tol2);
    mnt_vectorinterp_getVector(&vp, &data[0], &resVec[0]);
    std::cout << "at point " << target << " vector is " << resVec << '\n';

    // destroy
    mnt_vectorinterp_del(&vp);
    cloc->Delete();
    mnt_grid_del(&mgrid);
}

int main(int argc, char** argv) {

    testContainsPoint(10, 5);

    return 0;
}
