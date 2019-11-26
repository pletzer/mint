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
    // one or more ":" to discriminate file and mesh names
    ier = mnt_grid_loadFrom2DUgrid(&grd, "${CMAKE_SOURCE_DIR}/data/cs_16.nc::physics");
    assert(ier == 0);
    ier = mnt_grid_print(&grd);
    assert(ier == 0);
    size_t numCells;
    ier = mnt_grid_getNumberOfCells(&grd, &numCells);
    assert(ier == 0);
    std::cout << numCells << " cells\n";

    std::vector<vtkIdType> cellIds{0, 1, 2, 3, 12, 456, 1535};

    vtkIdType nodeIds[2];
    size_t edgeId;
    int signEdge;
    for (auto cellId : cellIds) {
        for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex) {
            ier = mnt_grid_getNodeIds(&grd, cellId, edgeIndex, nodeIds);
            std::cout << " cell " << cellId << " edge index " << edgeIndex 
                      << " connects to nodes "  << nodeIds[0] << ',' << nodeIds[1] << '\n';
            assert(ier == 0);
            ier = mnt_grid_getEdgeId(&grd, cellId, edgeIndex, &edgeId, &signEdge);
            std::cout << " cell " << cellId << " edge index " << edgeIndex 
                      << " maps to edge Id = " << edgeId << " (direction " << signEdge << ")\n";
            assert(ier == 0);
        }
    }

    // attach a field
    std::vector<double> cellByCellData(numCells * 4);
    for (size_t cellId = 0; cellId < numCells; ++cellId) {
        for (int ie = 0; ie < 4; ++ie) {
            ier = mnt_grid_getEdgeId(&grd, cellId, ie, &edgeId, &signEdge);
            assert(ier == 0);
            size_t k = cellId*4 + ie;
            cellByCellData[k] = (double) edgeId * (double) signEdge;        }
    }
    const char fieldName[] = "edge_ids";
    std::cout << "attaching edge field " << fieldName << " to grid\n";
    ier = mnt_grid_attach(&grd, fieldName, 4, &cellByCellData[0]);
    assert(ier == 0);
    
    ier = mnt_grid_dump(&grd, "cs_16.vtk");
    assert(ier == 0);

    ier = mnt_grid_del(&grd);
    assert(ier == 0); 
}


int main(int argc, char** argv) {

	testUgrid();
    testVTK();

    return 0;
}
