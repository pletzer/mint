#include <mntUgrid2D.h>
#undef NDEBUG // turn on asserts


void test1() {

    int ier;
    std::string file = "@CMAKE_SOURCE_DIR@/data/cs_4.nc";

    Ugrid2D ugr;
    ier = ugr.load(file, "physics");
    assert(ier == 0);

    size_t numPoints = ugr.getNumberOfPoints();
    size_t numEdges = ugr.getNumberOfEdges();
    size_t numFaces = ugr.getNumberOfFaces();

    assert(numFaces == 96);
    assert(numEdges == 192);
    assert(numPoints == 98);
}



int main() {

    test1();

    return 0;
}   
