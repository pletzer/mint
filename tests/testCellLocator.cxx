#include <mntCellLocator.h>
#include <cassert>
#undef NDEBUG // turn on asserts

void test1() {
    CellLocator_t* cloc;
    mnt_celllocator_new(&cloc);
    int num_cells = 1;
    const double verts[] = {0., 0., 0.,
                            1., 0., 0.,
                            1., 1., 0.,
                            0., 1., 0.};
    mnt_celllocator_setpoints(&cloc, num_cells, verts);
    mnt_celllocator_build(&cloc);
    const double target_point[] = {0.1, 0.2, 0.};
    double pcoords[3];
    long long cell_id;
    mnt_celllocator_find(&cloc, target_point, &cell_id, pcoords);
    std::cout << "cell id: " << cell_id << '\n';
    assert(cell_id == 0);
    mnt_celllocator_del(&cloc);
}


int main(int argc, char** argv) {

    test1();

    return 0;
}
