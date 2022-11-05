#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>

#include <mntGrid.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>



void saveEdgeVectors(Grid_t* grd, 
                 const std::vector<double>& u,
                 const std::vector<double>& v,
                 const std::string& filename) {
        
    std::size_t numEdges;
    mnt_grid_getNumberOfEdges(&grd, &numEdges);
    std::size_t numCells;
    mnt_grid_getNumberOfCells(&grd, &numCells);

    vtkUnstructuredGrid* vg = vtkUnstructuredGrid::New();
    vtkPoints* vp = vtkPoints::New();
    vtkDoubleArray* va = vtkDoubleArray::New();
    va->SetNumberOfComponents(3);
    va->SetNumberOfTuples(numEdges);
    vtkDoubleArray* vv = vtkDoubleArray::New();
    vv->SetName("vectors");
    vv->SetNumberOfComponents(3);
    vv->SetNumberOfTuples(numEdges);
    std::size_t edgeId;
    int edgeSign;
    Vec3 p0, p1, pe, pv;
    for (auto icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &edgeSign);
            mnt_grid_getPoints(&grd, icell, ie, &p0[0], &p1[0]);
            pe = 0.5*(p0 + p1);
            va->SetTuple(edgeId, &pe[0]);
            pv[0] = u[edgeId];
            pv[1] = v[edgeId];
            pv[2] = 0.0;
            vv->SetTuple(edgeId, &pv[0]);
        }
    }
    vp->SetData(va);
    vg->SetPoints(vp);
    vg->Allocate(numEdges);
    vtkIdList* ptIds = vtkIdList::New();
    ptIds->SetNumberOfIds(1);
    for (auto i = 0; i < numEdges; ++i) {
        ptIds->SetId(0, i);
        vg->InsertNextCell(VTK_VERTEX, ptIds);
    }
    ptIds->Delete();

    vg->GetPointData()->AddArray(vv);
    vtkUnstructuredGridWriter* vw = vtkUnstructuredGridWriter::New();
    vw->SetFileVersion(42);
    vw->SetInputData(vg);
    vw->SetFileName(filename.c_str());
    vw->Update();

    vw->Delete();
    vv->Delete();
    va->Delete();
    vp->Delete();
    vg->Delete();
}
