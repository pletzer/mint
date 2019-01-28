#include <mntGrid.h>
#include <vtkCellData.h>
#include <fstream>
#include <netcdf.h>
#include <string>
#include <cstring>

extern "C" 
int mnt_grid_new(Grid_t** self) {

    *self = new Grid_t();
    (*self)->verts = NULL;
    (*self)->pointData = NULL;
    (*self)->points = NULL;
    (*self)->grid = NULL;
    (*self)->reader = NULL;
    (*self)->writer = NULL;
    (*self)->doubleArrays.resize(0);
    return 0;
}

extern "C"
int mnt_grid_del(Grid_t** self) {

    for (size_t i = 0; i < (*self)->doubleArrays.size(); ++i) {
        (*self)->doubleArrays[i]->Delete();
    }
    if ((*self)->writer) (*self)->writer->Delete();
    if ((*self)->reader) {
        (*self)->reader->Delete();
    }
    else {
        if ((*self)->grid) (*self)->grid->Delete();
    }
    if ((*self)->points) (*self)->points->Delete();
    if ((*self)->pointData) (*self)->pointData->Delete();

    if ((*self)->verts) delete[] (*self)->verts;

    delete *self;

    return 0;
}

extern "C"
int mnt_grid_setPointsPtr(Grid_t** self, int nVertsPerCell, vtkIdType ncells, const double points[]) {

    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();

    int save = 1;
    int npoints = nVertsPerCell * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetNumberOfComponents(3);
    (*self)->pointData->SetVoidArray((double*) points, npoints*3, save);

    (*self)->points->SetData((*self)->pointData);

    (*self)->grid->Allocate(ncells, 1);

    int cellType = -1;
    if (nVertsPerCell == 4) {
        cellType = VTK_QUAD;
    }
    else if (nVertsPerCell == 8) {
        cellType = VTK_HEXAHEDRON;
    }
    else {
        // error
        return 1;
    }

    // connect
    vtkIdList* ptIds = vtkIdList::New();

    ptIds->SetNumberOfIds(nVertsPerCell);
    for (int i = 0; i < ncells; ++i) {
        for (int j = 0; j < nVertsPerCell; ++j) {
            ptIds->SetId(j, nVertsPerCell*i + j);
        }
        (*self)->grid->InsertNextCell(cellType, ptIds);
    }
    (*self)->grid->SetPoints((*self)->points);
    (*self)->grid->BuildLinks(); // DO WE NEED THIS?

    // clean
    ptIds->Delete();

    return 0;
}

extern "C"
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]) {

    if (!(*self)->grid) {
        return 1;
    }

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();

    vtkDoubleArray* vtkdata = vtkDoubleArray::New();
    vtkdata->SetName(varname);
    vtkdata->SetNumberOfTuples(ncells);
    vtkdata->SetNumberOfComponents(nDataPerCell);
    int save = 1;
    vtkdata->SetVoidArray((double*) data, ncells*nDataPerCell, save);

    // store
    (*self)->doubleArrays.push_back(vtkdata);

    // add to the grid
    (*self)->grid->GetCellData()->AddArray(vtkdata);

    return 0;
}


extern "C"
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr) {
    *grid_ptr = (*self)->grid;
    return 0;
}

extern "C"
int mnt_grid_loadFrom2DUgrid(Grid_t** self, const char* filename) {

    // open the file
    int ncid;
    int ier = nc_open(filename, NC_NOWRITE, &ncid);

    size_t ncells = 0;
    const int four = 4;

    size_t nlats = 0;
    size_t nlons = 0;
    double* lats = NULL;
    double* lons = NULL;
    int* quad_connectivity = NULL;

    // get the number of variables
    int nvars;
    ier = nc_inq_nvars(ncid, &nvars);

    // find the latitudes, longitudes and cell connectivity
    int latId = -1;
    int lonId = -1;
    char varname[NC_MAX_NAME];
    char standard_name[NC_MAX_NAME];
    char long_name[NC_MAX_NAME];
    char cf_role[NC_MAX_NAME];
    char units[NC_MAX_NAME];
    nc_type xtype;
    int ndims;
    int natts;
    int ntot;
    for (int ivar = 0; ivar < nvars; ++ivar) {

        // get the number of dimensions of this variable
        ier = nc_inq_varndims(ncid, ivar, &ndims);

        // get the dimensions of this variable
        int* dimids = new int[ndims];

        ier = nc_inq_var(ncid, ivar, varname, &xtype, &ndims, dimids, &natts);

        // get the variable attributes
        int ier1 = nc_get_att_text(ncid, ivar, "standard_name", standard_name);
        int ier2 = nc_get_att_text(ncid, ivar, "long_name", long_name);
        int ier3 = nc_get_att_text(ncid, ivar, "cf_role", cf_role);

        if (ier1 == NC_NOERR && ier2 == NC_NOERR) {
            // variable has "standard_name" and long_name attributes

            // is it a latitude or a longitude?
            if (strcmp(standard_name, "latitude") == 0 && 
                std::string(long_name).find("node") != std::string::npos) {

                // get the number of elements
                nlats = 1;
                for (int i = 0; i < ndims; ++i) {
                    size_t dim;
                    ier = nc_inq_dimlen(ncid, dimids[i], &dim);
                    nlats *= dim;
                }

                if (xtype == NC_DOUBLE) {

                    // allocate the data to receive the lats
                    lats = new double[ntot];

                    // read the latitudes
                    ier = nc_get_var_double(ncid, ivar, lats);
                }
            }
            else if (strcmp(standard_name, "longitude") == 0 && 
                std::string(long_name).find("node") != std::string::npos) {

                // get the number of elements
                nlons = 1;
                for (int i = 0; i < ndims; ++i) {
                    size_t dim;
                    ier = nc_inq_dimlen(ncid, dimids[i], &dim);
                    nlons *= dim;
                }

                if (xtype == NC_DOUBLE) {
                    // allocate the data to receive the lons
                    lons = new double[ntot];

                    // read the longitudes
                    ier = nc_get_var_double(ncid, ivar, lons);
                }
            }
        }
        else if (ier3 == NC_NOERR && 
                 strcmp(cf_role, "face_node_connectivity") == 0) {

            // get the number of cells
            ier = nc_inq_dimlen(ncid, dimids[0], &ncells);
            ntot = ncells * four;

            if (xtype == NC_INT) {

                // allocate the data
                quad_connectivity = new int[ntot];

                // read the connectivity
                ier = nc_get_var_int(ncid, ivar, quad_connectivity);
            }

        }

        delete[] dimids;
    }

    // close the netcdf file
    ier = nc_close(ncid);

    if (lons && lats && quad_connectivity) {

        // allocate the vertices and set the values
        (*self)->verts = new double[ncells * four * 3];
        for (int icell = 0; icell < ncells; ++icell) {
            for (int node = 0; node < four; ++node) {
                int k = quad_connectivity[node + four*icell];
                // even in 2d we have three components
                (*self)->verts[0 + icell*four*3] = lats[k];
                (*self)->verts[1 + icell*four*3] = lons[k];
                (*self)->verts[2 + icell*four*3] = 0.0;
            }
        }
    }


    // clean up
    if (lats) delete[] lats;
    if (lons) delete[] lons;
    if (quad_connectivity) delete[] quad_connectivity;

    return 0;
}

extern "C"
int mnt_grid_load(Grid_t** self, const char* filename) {
    // check if the file exists
    if (!fstream(filename).good()) {
        std::cerr << "ERROR file " << filename << " does not exist\n";
        return 1;        
    }

    if ((*self)->grid) {
        (*self)->grid->Delete();
    }
    (*self)->reader = vtkUnstructuredGridReader::New();
    (*self)->reader->SetFileName(filename);
    (*self)->reader->Update();
    (*self)->grid = (*self)->reader->GetOutput();
    return 0;
}

extern "C"
int mnt_grid_dump(Grid_t** self, const char* filename) {
    (*self)->writer = vtkUnstructuredGridWriter::New();
    (*self)->writer->SetFileName(filename);
    (*self)->writer->SetInputData((*self)->grid);
    (*self)->writer->Update();
    return 0;
}
