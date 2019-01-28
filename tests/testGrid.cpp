#include <mntGrid.h>
#undef NDEBUG // turn on asserts
#include <cassert>

void testVTK() {
    Grid_t* grd;
    mnt_grid_new(&grd);
    mnt_grid_load(&grd, "${CMAKE_SOURCE_DIR}/data/cs.vtk");
    mnt_grid_del(&grd);
}

void testUgrid() {
	int ier;
    Grid_t* grd;
    ier = mnt_grid_new(&grd);
    assert(ier == 0);
    ier = mnt_grid_loadFrom2DUgrid(&grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc");
    assert(ier == 0);
    ier = mnt_grid_del(&grd);
    assert(ier == 0); 
}


int main(int argc, char** argv) {

	testUgrid();
    testVTK();

    return 0;
}
