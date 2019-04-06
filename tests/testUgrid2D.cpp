#include <mntUgrid2D.h>
#undef NDEBUG // turn on asserts


void test1() {

    int ier;
    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ugr;
    ier = ugr.load(file, "physics");
    assert(ier == 0);

}



int main() {

    test1();

    return 0;
}   
