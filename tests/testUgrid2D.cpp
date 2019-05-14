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

void test2() {

  int ier;
  std::string file = "@CMAKE_SOURCE_DIR@/data/output_UGRID.nc";

  Ugrid2D ugr;
  ier = ugr.load(file, "Mesh2d");
  assert(ier == 0);
  ugr.dumpGridVtk("output.vtk");
}



int main() {

    test1();
    test2();

    return 0;
}   
