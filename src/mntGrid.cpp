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
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return 1;
    }

    size_t ncells = 0;
    const int four = 4;

    size_t nlats = 0;
    size_t nlons = 0;
    std::vector<double> lats;
    std::vector<double> lons;
    std::vector<vtkIdType> quad_connectivity;

    // get the number of variables
    int nvars;
    ier = nc_inq_nvars(ncid, &nvars);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: after inquiring the number of variables (ier = " 
                  << ier << ")\n";
        return 1;
    }

    // find the latitudes, longitudes and cell connectivity
    int latId = -1;
    int lonId = -1;
    char varname[NC_MAX_NAME];
    std::string standard_name(NC_MAX_NAME, ' ');
    std::string long_name(NC_MAX_NAME, ' ');
    std::string cf_role(NC_MAX_NAME, ' ');
    char units[NC_MAX_NAME];
    nc_type xtype;
    int ndims;
    int natts;
    int startIndex = 0;
    for (int ivar = 0; ivar < nvars; ++ivar) {

        // get the number of dimensions of this variable
        ier = nc_inq_varndims(ncid, ivar, &ndims);
        if (ier != NC_NOERR) {
            std::cerr << "ERROR: after inquiring the number of dimensions for var = " 
                      << ivar << " (ier = " << ier << ")\n";
            return 2;
        }

        // get the dimensions of this variable
        int* dimids = new int[ndims];

        ier = nc_inq_var(ncid, ivar, varname, &xtype, &ndims, dimids, &natts);

        // reset
        standard_name.assign(NC_MAX_NAME, ' ');
        long_name.assign(NC_MAX_NAME, ' ');
        cf_role.assign(NC_MAX_NAME, ' ');

        // get the variable attributes
        int ier1 = nc_get_att_text(ncid, ivar, "standard_name", &standard_name[0]);
        int ier2 = nc_get_att_text(ncid, ivar, "long_name", &long_name[0]);
        int ier3 = nc_get_att_text(ncid, ivar, "cf_role", &cf_role[0]);
        int ier4 = nc_get_att_int(ncid, ivar, "start_index", &startIndex);

        if (ier1 == NC_NOERR && ier2 == NC_NOERR) {

            // variable has "standard_name" and long_name attributes

            // get the number of elements
            int nelems = 1;
            for (int i = 0; i < ndims; ++i) {
                size_t dim;
                ier = nc_inq_dimlen(ncid, dimids[i], &dim);
                if (ier != NC_NOERR) {
                    std::cerr << "ERROR: after getting the dimension size (ier = " << ier << ")\n";
                        return 1;
                }
                    nelems *= dim;
            }

            // is it a latitude or a longitude?
            if (standard_name.find("latitude") != std::string::npos &&
                long_name.find("node") != std::string::npos) {

                // allocate the data to receive the lats
                lats.resize(nelems);

                if (xtype == NC_DOUBLE) {
                    // read the latitudes as doubles
                    ier = nc_get_var_double(ncid, ivar, &lats[0]);
                    if (ier != NC_NOERR) {
                        std::cerr << "ERROR: after reading latitudes as doubles (ier = " << ier << ")\n";
                        return 1;
                    }
                }
                else if (xtype == NC_FLOAT) {
                    // read the latitudes as floats
                    std::vector<float> data(nelems);
                    ier = nc_get_var_float(ncid, ivar, &data[0]);
                    if (ier != NC_NOERR) {
                        std::cerr << "ERROR: after reading latitudes as floats (ier = " << ier << ")\n";
                        return 1;
                    }
                    // copy 
                    for (size_t i = 0; i < nelems; ++i) {
                        lats[i] = (double) data[i];
                    }

                }
                else {

                }
            }
            else if (standard_name.find("longitude") != std::string::npos && 
                     long_name.find("node") != std::string::npos) {

                // allocate to receive the lons
                lons.resize(nelems);

                if (xtype == NC_DOUBLE) {
                    // read the longitudes as doubles
                    ier = nc_get_var_double(ncid, ivar, &lons[0]);
                    if (ier != NC_NOERR) {
                        std::cerr << "ERROR: after reading longitudes as double (ier = " << ier << ")\n";
                        return 1;
                    }
                }
                else if (xtype == NC_FLOAT) {
                    // read the longitudes as floats
                    std::vector<float> data(nelems);
                    ier = nc_get_var_float(ncid, ivar, &data[0]);
                    if (ier != NC_NOERR) {
                        std::cerr << "ERROR: after reading longitudes as floats (ier = " << ier << ")\n";
                        return 1;
                    }
                    // copy into double array
                    for (size_t i = 0; i < nelems; ++i) {
                        lons[i] = (double) data[i];
                    }
                }
            }
        }
        else if (ier3 == NC_NOERR && ier4 == NC_NOERR &&
                 cf_role.find("face_node_connectivity") != std::string::npos) {

            // get the number of cells
            ier = nc_inq_dimlen(ncid, dimids[0], &ncells);
            size_t nelems = ncells * four;

            // allocate the data
            quad_connectivity.resize(nelems);

            // read the connectivity
            if (xtype == NC_INT) {

                std::vector<int> data(nelems);
                ier = nc_get_var_int(ncid, ivar, &data[0]);
                if (ier != NC_NOERR) {
                    std::cerr << "ERROR: after reading cell connectivity using ints (ier = " << ier << ")\n";
                    return 1;
                }
                // copy into vtkIdType array
                for (size_t i = 0; i < nelems; ++i) {
                    quad_connectivity[i] = (vtkIdType) data[i] - startIndex;
                }

            }
            else if (xtype == NC_INT64) {

                ier = nc_get_var_longlong(ncid, ivar, &quad_connectivity[0]);
                if (ier != NC_NOERR) {
                    std::cerr << "ERROR: after reading cell connectivity using int64 (ier = " << ier << ")\n";
                    return 1;
                }
                // substract start index
                for (size_t i = 0; i < nelems; ++i) {
                    quad_connectivity[i] -= startIndex;
                }
            }

        }

        delete[] dimids;
    }

    // close the netcdf file
    ier = nc_close(ncid);

    if (lons.size() > 0 && lats .size() > 0 && quad_connectivity.size() > 0) {

        // allocate the vertices and set the values
        (*self)->verts = new double[ncells * four * 3];

        std::vector<double> diffLonMinusZeroPlus(four);

        for (size_t icell = 0; icell < ncells; ++icell) {

            // fix longitude if crossing the dateline
            // use the first longitude as the base
            size_t kBase = quad_connectivity[four*icell];
            double lonBase = lons[kBase];


            for (int node = 0; node < four; ++node) {

                size_t k = quad_connectivity[node + four*icell];
                double lon = lons[k];

                // add/subtract 360.0, whatever it takes to reduce the distance 
                // between this longitude and the base longitude
                double diffLon = lon - lonBase;
                diffLonMinusZeroPlus[0] = std::abs(diffLon - 360.);
                diffLonMinusZeroPlus[1] = std::abs(diffLon - 0.);
                diffLonMinusZeroPlus[2] = std::abs(diffLon + 360.);
                std::vector<double>::iterator it = std::min_element(diffLonMinusZeroPlus.begin(), 
                                                                    diffLonMinusZeroPlus.end());
                int indexMin = (int) std::distance(diffLonMinusZeroPlus.begin(), it);
                // fix the longitude
                lon += (indexMin - 1) * 360.0;

                // even in 2d we have three components
                (*self)->verts[0 + node*3 + icell*four*3] = lats[k];
                (*self)->verts[1 + node*3 + icell*four*3] = lon;
                (*self)->verts[2 + node*3 + icell*four*3] = 0.0;
            }
        }
    }

    // set the pointer
    ier = mnt_grid_setPointsPtr(self, (vtkIdType) four, (vtkIdType) ncells, (*self)->verts);

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

extern "C"
int mnt_grid_print(Grid_t** self) {

    vtkPoints* points = (*self)->grid->GetPoints();
    vtkIdType npoints = points->GetNumberOfPoints();
    std::cerr << "Number of points: " << npoints << '\n';

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();
    std::cerr << "Number of cells: " << ncells << '\n';

    std::vector<double> pt(3);

    for (vtkIdType i = 0; i < ncells; ++i) {

        vtkCell* cell = (*self)->grid->GetCell(i);

        for (int j = 0; j < cell->GetNumberOfPoints(); ++j) {
            vtkIdType k = cell->GetPointId(j);
            (*self)->points->GetPoint(k, &pt[0]);
            std::cout << "\tpoint " << pt[0] << ',' << pt[1] << ',' << pt[2] << '\n';
        }
    }

    return 0;
}
