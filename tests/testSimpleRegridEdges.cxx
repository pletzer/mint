#include "mntRegridEdges.h"
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIdList.h>
#include <iostream>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>



int main() {

    double v0[] = {0., 0., 0.};
    double v1[] = {1., 0., 0.};
    double v2[] = {1., 1., 0.};
    double v3[] = {0., 1., 0.};

    vtkIdList* ptIds;

    // construct the src grid
    vtkPoints* srcPoints = vtkPoints::New();
    srcPoints->SetDataTypeToDouble();
    srcPoints->InsertNextPoint(v0);
    srcPoints->InsertNextPoint(v1);
    srcPoints->InsertNextPoint(v2);
    srcPoints->InsertNextPoint(v3);
    vtkUnstructuredGrid* srcGrid = vtkUnstructuredGrid::New();
    srcGrid->SetPoints(srcPoints);
    srcGrid->Allocate(1, 1);
    ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    srcGrid->InsertNextCell(VTK_QUAD, ptIds);
    ptIds->Delete();

    // construct dst grid
    vtkPoints* dstPoints = vtkPoints::New();
    dstPoints->SetDataTypeToDouble();
    dstPoints->InsertNextPoint(v0);
    dstPoints->InsertNextPoint(v1);
    dstPoints->InsertNextPoint(v2);
    dstPoints->InsertNextPoint(v3);
    vtkUnstructuredGrid* dstGrid = vtkUnstructuredGrid::New();
    dstGrid->SetPoints(dstPoints);
    dstGrid->Allocate(1, 1);
    ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(4);
    for (vtkIdType i = 0; i < 4; ++i) {
        ptIds->SetId(i, i);
    }
    dstGrid->InsertNextCell(VTK_QUAD, ptIds);
    ptIds->Delete();

    RegridEdges_t* rg;
    int ier;

    ier = mnt_regridedges_new(&rg);
    assert(ier == 0);

    ier = mnt_regridedges_setSrcGrid(&rg, srcGrid);
    assert(ier == 0);

    ier = mnt_regridedges_setDstGrid(&rg, dstGrid);
    assert(ier == 0);

    int num_cells_per_bucket = 8;
    ier = mnt_regridedges_build(&rg, num_cells_per_bucket);
    assert(ier == 0);
    
    std::string outputFilename = "simpleRegridEdgesWeights.nc";
    ier = mnt_regridedges_dump(&rg, outputFilename.c_str(), outputFilename.size());
    assert(ier == 0);

    double srcData[] = {1., 2., 1., 2.};
    double dstData[4];

    ier = mnt_regridedges_applyWeights(&rg, srcData, dstData);
    assert(ier == 0);

    ier = mnt_regridedges_del(&rg);
    assert(ier == 0);

    for (size_t i = 0; i < 4; ++i) {
        std::cout << "edge " << i << " src, dst value = " << srcData[i] << ',' << dstData[i] << '\n';
    }

    // clean up
    srcGrid->Delete();

}