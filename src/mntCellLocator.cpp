#include <mntCellLocator.h>
#include <limits>
#include <cstring>
#include <string>

extern "C"
int mnt_celllocator_new(CellLocator_t** self) {
    *self = new CellLocator_t();
    mnt_grid_new(&(*self)->gridt);
    (*self)->loc = vtkCellLocator::New();
    (*self)->cell = vtkGenericCell::New();
    return 0;
}

extern "C"
int mnt_celllocator_load(CellLocator_t** self, const char* fort_filename, size_t n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    // load the grid
    int ier = mnt_grid_load(&(*self)->gridt, filename.c_str());
    if (ier != 0) {
        return ier;
    }
    return 0;
}

extern "C"
int mnt_celllocator_del(CellLocator_t** self) {
    (*self)->cell->Delete();
    (*self)->loc->Delete();
    mnt_grid_del(&(*self)->gridt);
    delete *self;
    return 0;
}

extern "C"
int mnt_celllocator_setPointsPtr(CellLocator_t** self, int nVertsPerCell, size_t ncells, 
		                         const double points[]) {
    int ier = mnt_grid_setPointsPtr(&(*self)->gridt, nVertsPerCell, ncells, points);
    return ier;
}

extern "C"
int mnt_celllocator_build(CellLocator_t** self, int num_cells_per_bucket) {

    // get the grid
    vtkUnstructuredGrid* grd;
    mnt_grid_get(&(*self)->gridt, &grd);

    // build the locator
    (*self)->loc->SetNumberOfCellsPerBucket(num_cells_per_bucket);
    (*self)->loc->SetDataSet(grd);
    (*self)->loc->BuildLocator();

    return 0;
}

extern "C"
int mnt_celllocator_checkGrid(CellLocator_t** self, double tol, int* numBadCells) {

    // initialize
    *numBadCells = 0;
    int ier = 0;

    vtkIdType ncells = (*self)->gridt->grid->GetNumberOfCells();
    if (ncells == 0) {
        // no op if there are no cells
        return 0;
    }

    // local cell indices of vertices used to compute the 4 tet volumes
    const vtkIdType j0[] = {1, 3, 4, 6};
    const vtkIdType ja[] = {2, 0, 7, 5};
    const vtkIdType jb[] = {0, 2, 5, 7};
    vtkIdType jc[] = {5, 7, 0, 2};

    // cell vertices
    vtkPoints* points = (*self)->gridt->grid->GetPoints();
    // four vertices spanning a tetrahedron
    double p0[3], pa[3], pb[3], pc[3];

    // vectors from the base vertex to other 3 verts
    double a[3], b[3], c[3];

    size_t numElems = 4; // number of tets
    double isQuad = 0.;

    unsigned char cellType = (*self)->gridt->grid->GetCell(0)->GetCellType();
    if (cellType == VTK_QUAD) {
        // fake vertical dimension
        c[0] = 0.;
        c[1] = 0.;
        c[2] = 1.;
        jc[0] = 0; jc[1] = 0; // won't use those must be sure to create mem issues
        numElems = 2; // number of triangles
        isQuad = 1.;
    } 
    else if (cellType != VTK_HEXAHEDRON) {
        // not supported 
        return 1;
    } 

    // iterate over the cells
    for (vtkIdType icell = 0; icell < ncells; ++icell) {

        // number of negative area/volumes in each cell
        int numBad = 0;
        // vertex Ids for that cell
        vtkIdList* ptIds = (*self)->gridt->grid->GetCell(icell)->GetPointIds();

        // iterate over the 2 triangles or 4 tets
        for (size_t j = 0; j < numElems; ++j) {

            // get the point Ids for each vertex of the triangle/tet
            vtkIdType k0 = ptIds->GetId(j0[j]);
            vtkIdType ka = ptIds->GetId(ja[j]);
            vtkIdType kb = ptIds->GetId(jb[j]);
            vtkIdType kc = ptIds->GetId(jc[j]);
            points->GetPoint(k0, p0);
            points->GetPoint(ka, pa);
            points->GetPoint(kb, pb);
            points->GetPoint(kc, pc);

            // vectors from the base point to the other triangle/tet vertices
            for (size_t dim = 0; dim < 3; ++dim) {
                a[dim] = pa[dim] - p0[dim];
                b[dim] = pb[dim] - p0[dim];
                // use the fake vertical (0, 0, 1) vector for quads
                c[dim] = isQuad*c[dim] + (1. - isQuad)*(pc[dim] - p0[dim]);
            }

            // area or volume
            double volume = (a[1]*b[2] - a[2]*b[1])*c[0] 
                          + (a[2]*b[0] - a[0]*b[2])*c[1]
                          + (a[0]*b[1] - a[1]*b[0])*c[2];

            if (volume < tol) {
                numBad++;
                std::cout << "Bad cell (0-based Id = " << icell << ") point Ids: " 
                          << j0[j] << ',' << ja[j] << ',' << jb[j] << ',' << jc[j] << " coords: " 
                          << p0[0] << ',' << p0[1] << ',' << p0[2] << ';' 
                          << pa[0] << ',' << pa[1] << ',' << pa[2] << ';' 
                          << pb[0] << ',' << pb[1] << ',' << pb[2] << ';' 
                          << pc[0] << ',' << pc[1] << ',' << pc[2] << '\n';
            }
        }

        if (numBad > 0) {
            (*numBadCells)++;
        }
    }

    return ier;
}

extern "C"
int mnt_celllocator_runGridDiagnostics(CellLocator_t** self) {
    vtkIdType ncells = (*self)->gridt->grid->GetNumberOfCells();
    vtkIdType npts = (*self)->gridt->grid->GetNumberOfPoints();
    std::cout << "No of cells = " << ncells << " no of points = " << npts << '\n';
    vtkPoints* points = (*self)->gridt->grid->GetPoints();
    vtkCell* cell;
    double vert[3];
    for (vtkIdType i = 0; i < ncells; ++i) {
        std::cout << "cell " << i << " verts ";
        cell = (*self)->gridt->grid->GetCell(i);
        vtkIdList* ptIds = cell->GetPointIds();
        for (vtkIdType j = 0; j < ptIds->GetNumberOfIds(); ++j) {
            vtkIdType k = ptIds->GetId(j);
            points->GetPoint(k, vert);
            std::cout << "\t" << k << " (";
            for (int el = 0; el < 3; ++el) {
                std::cout << vert[el] << ',';
            }
            std::cout << ") ";
        }
        std::cout << '\n';
    }
    return 0;
}


extern "C"
int mnt_celllocator_dumpGrid(CellLocator_t** self, const char* fort_filename, size_t n) {
    // copy fortran string into c string
    std::string filename = std::string(fort_filename, n);
    int ier = mnt_grid_dump(&(*self)->gridt, filename.c_str());
    return ier;
}

extern "C"
int mnt_celllocator_find(CellLocator_t** self, const double point[], long long* cellId, 
		                 double pcoords[]) {
    // may need to be adjusted
    const double tol2 = 10 * std::numeric_limits<double>::epsilon(); 
    *cellId = (*self)->loc->FindCell((double*) point, tol2, (*self)->cell, pcoords, (*self)->weights);
    return 0;
}

extern "C"
int mnt_celllocator_interpPoint(CellLocator_t** self, long long cellId, const double pcoords[], double point[]) {
    if (cellId < 0) {
        return 1;
    }
    int subId = 0;
    vtkCell* cell = (*self)->gridt->grid->GetCell(cellId);
    cell->EvaluateLocation(subId, (double*) pcoords, point, (*self)->weights);
    return 0;
}

extern "C"
void mnt_celllocator_printAddress(void* something) {
    printf("address is %lld (0x%16zx)\n", (long long) something, (size_t) something);
}
