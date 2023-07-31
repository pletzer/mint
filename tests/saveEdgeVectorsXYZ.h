#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>

#include <mntGrid.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>



void saveEdgeVectorsXYZ(Grid_t* grd, 
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
    Vec3 p0, p1, pe, pv, xyz;
    for (auto icell = 0; icell < numCells; ++icell) {
        for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {
            mnt_grid_getEdgeId(&grd, icell, ie, &edgeId, &edgeSign);
            mnt_grid_getPoints(&grd, icell, ie, &p0[0], &p1[0]);
            pe = 0.5*(p0 + p1);
            va->SetTuple(edgeId, &pe[0]);
            double lam = pe[0] * M_PI/180.;
            double the = pe[1] * M_PI/180.;
            double rho = cos(the);
            double x = rho*cos(lam);
            double y = rho*sin(lam);
            double z = sin(the);
            xyz[0] = x;
            xyz[1] = y;
            xyz[2] = z;
            va->SetTuple(edgeId, &xyz[0]);
            pv[0] = - u[edgeId]*sin(lam) - v[edgeId]*z*cos(lam);
            pv[1] = + u[edgeId]*cos(lam) - v[edgeId]*z*sin(lam);
            pv[2] = v[edgeId]*rho;
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
