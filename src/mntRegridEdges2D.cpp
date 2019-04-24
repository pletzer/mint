#include <mntRegridEdges2D.h>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <netcdf.h>

/**
 * Extract the file name and mesh name from non-zero terminated fortran string 
 * @param fort_filename eg "fileName:meshName"
 * @param n length o filename
 * @param fileName file name (output)
 * @param meshName mesh name (output)
 */
void getFileAndMeshNames(const char* fort_filename, int n, 
                        std::string& fileName, std::string& meshName) {

    std::string fileAndMeshName = std::string(fort_filename, n);

    size_t column = fileAndMeshName.find(':');
    fileName = fileAndMeshName.substr(0, column);

    if (column != std::string::npos) {
        meshName = fileAndMeshName.substr(column + 1);
    }
}

extern "C"
int mnt_regridedges2d_new(RegridEdges2D_t** self) {

    *self = new RegridEdges2D_t();
    (*self)->numSrcEdges = 0;
    (*self)->numDstEdges = 0;

    return 0;
}

extern "C"
int mnt_regridedges2d_del(RegridEdges2D_t** self) {

    delete *self;

    return 0;
}

extern "C"
int mnt_regridedges2d_loadEdgeField(RegridEdges2D_t** self,
                                    const char* fort_filename, int nFilenameLength,
                                    const char* field_name, int nFieldNameLength,
                                    size_t ndata, double data[]) {

    std::string fileAndMeshName = std::string(fort_filename, nFilenameLength);

    // filter out the mesh name, if present (not used here)
    size_t columnL = fileAndMeshName.find(':');
    std::string filename = fileAndMeshName.substr(0, columnL);

    std::string fieldname = std::string(field_name, nFieldNameLength);

    // open the file
    int ncid;
    int ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open \"" << filename << "\"\n";
        nc_close(ncid);
        return 1;
    }

    // check if the variable/field exists
    int varId;
    ier = nc_inq_varid(ncid, fieldname.c_str(), &varId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not find variable \"" << fieldname << "\"\n";
        nc_close(ncid);
        return 1;
    }

    // check that the field has the "location" attribute
    size_t nLoc;
    ier = nc_inq_attlen(ncid, varId, "location", &nLoc);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: variable \"" << fieldname << "\" does not appear to have attribute 'location' (ier = " << ier << ")\n";
        nc_close(ncid);
        return 2;
    }
    char location[nLoc + 1];
    ier = nc_get_att_text(ncid, varId, "location", location);
    location[nLoc] = '\0';
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: attribute \"location\" of variable \"" << fieldname << "\" could not be read (ier = " << ier << ")\n";
        nc_close(ncid);
        return 6;
    }
    // check location is set to "edge"
    if (strcmp(location, "edge") != 0) {
        std::cerr << "ERROR: attribute \"location\" of variable " << fieldname << " is not edge  ("
                  << location << ")\n";
        nc_close(ncid);
        return 3;
    }

    // check if the data has the right dimension
    int ndims;
    ier = nc_inq_varndims(ncid, varId, &ndims);
    int dimIds[ndims];
    ier = nc_inq_vardimid(ncid, varId, dimIds);
    size_t n;
    ier = nc_inq_dimlen(ncid, dimIds[0], &n);
    if (n != ndata) {
        std::cerr << "ERROR: size of \"" << fieldname << "\" should be " << n
                  << " but got " << ndata << "\n";
        nc_close(ncid);
        return 5;        
    }

    // TO DO 
    // is there are way to check if a field is an edge integral vs a vector field? 
    // Assume field is a line integral

    // now read
    ier = nc_get_var_double(ncid, varId, data);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: while reading variable '" << fieldname << "'\n";
        nc_close(ncid);
        return 4;
    }

    // close the netcdf file
    ier = nc_close(ncid);

    return 0;
}

extern "C"
int mnt_regridedges2d_dumpEdgeField(RegridEdges2D_t** self,
                                    const char* fort_filename, int nFilenameLength,
                                    const char* field_name, int nFieldNameLength,
                                    size_t ndata, const double data[]) {
    
    std::string filename = std::string(fort_filename, nFilenameLength);
    std::string fieldname = std::string(field_name, nFieldNameLength);

    int ncid, ier;
    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not create file \"" << filename << "\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // create dimensions
    int numEdgesId;

    ier = nc_def_dim(ncid, "num_edges", ndata, &numEdgesId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_edges\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }    

    // create variable
    int dataId;
    int dims[] = {numEdgesId};
    ier = nc_def_var(ncid, fieldname.c_str(), NC_DOUBLE, 1, dims, &dataId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"data\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    // write the data
    ier = nc_put_var_double(ncid, dataId, data);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"data\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    // close the netcdf file
    ier = nc_close(ncid);    

    return 0;
}


extern "C"
int mnt_regridedges2d_loadSrcGrid(RegridEdges2D_t** self, 
		                          const char* fort_filename, int n) {

    int ier;
    std::string fileName, meshName;

    getFileAndMeshNames(fort_filename, n, fileName, meshName);

    if(fileName.size() == 0) return 1; // could not extract the file name
    if(meshName.size() == 0) return 2; // could not extract the mesh name

    ier = (*self)->srcGrid.load(fileName, meshName);

    (*self)->numSrcEdges = (*self)->srcGrid.getNumberOfEdges();

    return ier;
}

extern "C"
int mnt_regridedges2d_loadDstGrid(RegridEdges2D_t** self, 
		                          const char* fort_filename, int n) {

    int ier;
    std::string fileName, meshName;

    getFileAndMeshNames(fort_filename, n, fileName, meshName);

    if(fileName.size() == 0) return 1; // could not extract the file name
    if(meshName.size() == 0) return 2; // could not extract the mesh name

    ier = (*self)->dstGrid.load(fileName, meshName);

    (*self)->numDstEdges = (*self)->dstGrid.getNumberOfEdges();

    return ier;
}

extern "C"
int mnt_regridedges2d_build(RegridEdges2D_t** self, int numCellsPerBucket) {

    // checks
    if ((*self)->numSrcEdges == 0) {
        std::cerr << "mnt_regridedges_build: ERROR must load source grid!\n";
        return 1;
    }
    if ((*self)->numDstEdges == 0) {
        std::cerr << "mnt_regridedges_build: ERROR must load destination grid!\n";
        return 2;
    }

    (*self)->weightDstEdgeIds.resize(0);
    (*self)->weightSrcEdgeIds.resize(0);
    (*self)->weights.resize(0);
    // MAY WANT TO RESERVE SPACE

    Vector<double>  dstXi0(3, 0.0);
    Vector<double>  dstXi1(3, 0.0);
    Vector<double>  srcXi0(3, 0.0);
    Vector<double>  srcXi1(3, 0.0);

    // build the source grid locator
    (*self)->srcGrid.buildLocator(numCellsPerBucket);

    // compute the weights
    for (size_t dstEdgeId = 0; dstEdgeId < (*self)->numDstEdges; ++dstEdgeId) {

        // get the start end points of the dst egde
        std::vector< Vector<double> > dstEdgePoints = (*self)->dstGrid.getEdgePointsRegularized(dstEdgeId);
        Vector<double> u = dstEdgePoints[1] - dstEdgePoints[0];

        // find all the intersections between the dst edge and the source grid
        std::vector< std::pair<size_t, std::vector<double> > > intersections = 
                (*self)->srcGrid.findIntersectionsWithLine(dstEdgePoints[0], dstEdgePoints[1]);

        for (const std::pair<size_t, std::vector<double> >& srcCellIdLambdas : intersections) {

            size_t srcCellId = srcCellIdLambdas.first;

            // compute and set the regularized nodes of the cell
            (*self)->srcGrid.setCellPoints(srcCellId);

            double lambda0 = srcCellIdLambdas.second[0];
            double lambda1 = srcCellIdLambdas.second[1];
            Vector<double> dstPoint0 = dstEdgePoints[0] + lambda0*u;
            Vector<double> dstPoint1 = dstEdgePoints[0] + lambda1*u;

            // compute the src cell parametric coords of the dst edge segment
            bool inside;
            inside = (*self)->srcGrid.getParamCoords(dstPoint0, &dstXi0[0]); // need to check
            inside = (*self)->srcGrid.getParamCoords(dstPoint1, &dstXi1[0]); // need to check

            // iterate over the edges of the src cell
            const size_t* srcEdgeIds = (*self)->srcGrid.getFaceEdgeIds(srcCellId);
            for (size_t i = 0; i < 4; ++i) { // 2d (4 edges)
                size_t srcEdgeId = srcEdgeIds[i];
                // get the the end points of this src cell edge
                std::vector< Vector<double> > srcPoints = (*self)->srcGrid.getEdgePointsRegularized(srcEdgeId);
                const Vector<double>& srcPoint0 = srcPoints[0];
                const Vector<double>& srcPoint1 = srcPoints[1];
                // compute the src cell parametric coords of the src edge
                inside = (*self)->srcGrid.getParamCoords(srcPoint0, &srcXi0[0]); // need to check
                inside = (*self)->srcGrid.getParamCoords(srcPoint1, &srcXi1[0]); // need to check

                // compute the interpolation weight
                double weight = 1.0;
                for (size_t d = 0; d < 2; ++d) {

                    // mid point of target
                    double dstXiMid = 0.5*(dstXi0[d] + dstXi1[d]);
                    // length of target
                    double dstXiLen = dstXi1[d] - dstXi0[d];

                    // mid point of src edge in parameter space
                    double x = 0.5*(srcXi0[d] + srcXi0[d]);

                    // use Lagrange interpolation to evaluate the basis function integral for
                    // any for the 3 possible x values in {0, 0.5, 1}. This formula will make 
                    // it easier to extend the code to 3d
                    double xm00 = x;
                    double xm05 = x - 0.5;
                    double xm10 = x - 1.0;
                    double lag00 = + 2 * xm05 * xm10;
                    double lag05 = - 4 * xm00 * xm10;
                    double lag10 = + 2 * xm00 * xm05;

                    // CHEK FOR SIGN!!!
                    weight *= (1.0 - dstXiMid)*lag00 + dstXiLen*lag05 + dstXiMid*lag10;
                }

                (*self)->weightDstEdgeIds.push_back((long long) dstEdgeId);
                (*self)->weightSrcEdgeIds.push_back((long long) srcEdgeId);
                (*self)->weights.push_back(weight);
            }

        }

    }

    return 0;
}

extern "C"
int mnt_regridedges2d_getNumSrcEdges(RegridEdges2D_t** self, size_t* nPtr) {
    *nPtr = (*self)->numSrcEdges;
    return 0;
}

extern "C"
int mnt_regridedges2d_getNumDstEdges(RegridEdges2D_t** self, size_t* nPtr) {
    *nPtr = (*self)->numDstEdges;
    return 0;
}

extern "C"
int mnt_regridedges2d_apply(RegridEdges2D_t** self, 
	                        const double src_data[], double dst_data[]) {


    if ((*self)->numSrcEdges == 0 || (*self)->numDstEdges == 0) {
        std::cerr << "ERROR: looks like the src grid connectivity is not set.\n";
        std::cerr << "Typically this would occur if you did not read the grid\n";
        std::cerr << "from the netcdf Ugrid file.\n";
        return 1;
    }

    // initialize the dst data to zero
    for (size_t i = 0; i < (*self)->numDstEdges; ++i) {
        dst_data[i] = 0.0;
    }

    // add the contributions from each src edge
    for (size_t i = 0; i < (*self)->weights.size(); ++i) {

        size_t dstEdgeId = (size_t) (*self)->weightDstEdgeIds[i];
        size_t srcEdgeId = (size_t) (*self)->weightSrcEdgeIds[i];
        double weight = (*self)->weights[i];

        dst_data[dstEdgeId] += (*self)->weights[i] * src_data[srcEdgeId];

    }

    return 0;
}


extern "C"
int mnt_regridedges2d_loadWeights(RegridEdges2D_t** self, 
                                const char* fort_filename, int n) {

    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);

    int ncid, ier;
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not open file \"" << filename << "\"!\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // get the sizes
    size_t numWeights;
    int numWeightsId;
    ier = nc_inq_dimid(ncid, "num_weights", &numWeightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not inquire dimension \"num_weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }
    ier = nc_inq_dimlen(ncid, numWeightsId, &numWeights);


    int dstEdgeIdsId, srcEdgeIdsId, weightsId;

    ier = nc_inq_varid(ncid, "dst_edge_ids", &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }
    ier = nc_inq_varid(ncid, "src_edge_ids", &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }
    ier = nc_inq_varid(ncid, "weights", &weightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 7;
    }

    (*self)->weights.resize(numWeights);
    (*self)->weightDstEdgeIds.resize(numWeights);
    (*self)->weightSrcEdgeIds.resize(numWeights);

    // read
    ier = nc_get_var_double(ncid, weightsId, &((*self)->weights)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read var \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 8;
    }
    ier = nc_get_var_longlong(ncid, dstEdgeIdsId, &((*self)->weightDstEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_get_var_longlong(ncid, srcEdgeIdsId, &((*self)->weightSrcEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }

    ier = nc_close(ncid);    
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not close file \"" << filename << "\"!\n";
        std::cerr << nc_strerror (ier);
        return 13;
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_dumpWeights(RegridEdges2D_t** self, 
		                        const char* fort_filename, int n) {

    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);

    size_t numWeights = (*self)->weights.size();

    int ncid, ier;
    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not create file \"" << filename << "\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        return 1;
    }

    // create dimensions
    int numWeightsId;
    ier = nc_def_dim(ncid, "num_weights", (int) numWeights, &numWeightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_weights\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }

    // create variables
    int numWeightsAxis[] = {numWeightsId};

    int dstEdgeIdsId;
    ier = nc_def_var(ncid, "dst_edge_ids", NC_INT64, 1, numWeightsAxis, &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_edge_ids\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    int srcEdgeIdsId;
    ier = nc_def_var(ncid, "src_edge_ids", NC_INT64, 1, numWeightsAxis, &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    int weightsId;
    ier = nc_def_var(ncid, "weights", NC_DOUBLE, 1, numWeightsAxis, &weightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 7;
    }

    // close define mode
    ier = nc_enddef(ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not end define mode\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 8;
    }

    ier = nc_put_var_longlong(ncid, dstEdgeIdsId, &((*self)->weightDstEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"dst_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_put_var_longlong(ncid, srcEdgeIdsId, &((*self)->weightSrcEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"src_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_put_var_double(ncid, weightsId, &((*self)->weights)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"weights\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 12;
    }

    ier = nc_close(ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not close file \"" << filename << "\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 13;
    }

    return 0;
}

extern "C"
int mnt_regridedges2d_print(RegridEdges2D_t** self) {

    size_t numWeights = (*self)->weights.size();
    std::cout << "Number of weights: " << numWeights << '\n';

    printf("     index  dstEdgeId  srcEdgeId          weight\n");
    for (size_t i = 0; i < numWeights; ++i) {
        printf("%10ld %10ld %10ld %15.8le\n", i, (*self)->weightDstEdgeIds[i], 
            (*self)->weightSrcEdgeIds[i], (*self)->weights[i]);
    }

    return 0;
}

