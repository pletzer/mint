#include <mntGrid.h>

void testVTK() {
    Grid_t* grd;
    mnt_grid_new(&grd);
    mnt_grid_load(&grd, "${CMAKE_SOURCE_DIR}/data/cs.vtk");
    mnt_grid_del(&grd);
}

void testUgrid() {
    Grid_t* grd;
    mnt_grid_new(&grd);
    mnt_grid_loadFrom2DUgrid(&grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc");
    mnt_grid_del(&grd);
}


int main(int argc, char** argv) {

	testUgrid();
    testVTK();

    return 0;
}
