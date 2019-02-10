#include <mntRegridEdges.h>
#include <mntPolysegmentIter.h>
#include <mntQuadEdgeIter.h>
#include <iostream>
#include <cstdio>
#include <vtkIdList.h>
#include <netcdf.h>
#include <vtkHexahedron.h> // for 3d grids
#include <vtkQuad.h>       // for 2d grids
#include <vtkCell.h>

extern "C"
int mnt_regridedges_new(RegridEdges_t** self) {
    *self = new RegridEdges_t();
    (*self)->srcGrid = NULL;
    (*self)->dstGrid = NULL;
    (*self)->srcLoc = vtkCellLocator::New();
    (*self)->numSrcCells = 0;
    (*self)->numDstCells = 0;
    (*self)->numPointsPerCell = 4; // 2d
    (*self)->numEdgesPerCell = 4;  // 2d
    (*self)->srcGridObj = NULL;
    (*self)->dstGridObj = NULL;


    return 0;
}

extern "C"
int mnt_regridedges_del(RegridEdges_t** self) {

    // destroy the cell locator
    (*self)->srcLoc->Delete();
   
    // destroy the source and destination grids if this instance owns them
    if ((*self)->srcGridObj) {
        mnt_grid_del(&((*self)->srcGridObj));
    }
    if ((*self)->dstGridObj) {
        mnt_grid_del(&((*self)->dstGridObj));
    }

    delete *self;

    return 0;
}

extern "C"
int mnt_regridedges_loadSrc(RegridEdges_t** self, 
		            const char* fort_filename, int n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    if (!(*self)->srcGridObj) {
        mnt_grid_new(&((*self)->srcGridObj));
    }
    int ier = mnt_grid_loadFrom2DUgrid(&(*self)->srcGridObj, filename.c_str());
    (*self)->srcGrid = (*self)->srcGridObj->grid;
    return ier;
}

extern "C"
int mnt_regridedges_loadDst(RegridEdges_t** self, 
		            const char* fort_filename, int n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    if (!(*self)->dstGridObj) {
        mnt_grid_new(&((*self)->dstGridObj));
    }
    int ier = mnt_grid_loadFrom2DUgrid(&(*self)->dstGridObj, filename.c_str());
    (*self)->dstGrid = (*self)->dstGridObj->grid;
    return ier;
}


extern "C"
int mnt_regridedges_setSrcGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid) {
    // check
    if (grid->GetNumberOfPoints() != (*self)->numPointsPerCell * grid->GetNumberOfCells()) {
        std::cerr << "mnt_regridedges_setSrcGrid: ERROR num points != "
                  << (*self)->numPointsPerCell << " * num cells!\n";
        return 1;
    }

    // borrow pointer
    (*self)->srcGrid = grid;

    return 0;
}

extern "C"
int mnt_regridedges_setDstGrid(RegridEdges_t** self, vtkUnstructuredGrid* grid) {
    // check
    if (grid->GetNumberOfPoints() != (*self)->numPointsPerCell * grid->GetNumberOfCells()) {
        std::cerr << "mnt_regridedges_setDstGrid: ERROR num points != "
                  << (*self)->numPointsPerCell << "* num cells!\n";
        return 1;
    }

    // borrow pointer
    (*self)->dstGrid = grid;
    return 0;
}

extern "C"
int mnt_regridedges_setSrcPointsPtr(RegridEdges_t** self, size_t nVertsPerCell, size_t ncells, const double verts[]) {
    mnt_grid_new(&(*self)->srcGridObj);
    mnt_grid_setPointsPtr(&(*self)->srcGridObj, (int) nVertsPerCell, (vtkIdType) ncells, verts);
    mnt_grid_get(&(*self)->srcGridObj, &(*self)->srcGrid);
    return 0;
}

extern "C"
int mnt_regridedges_setDstPointsPtr(RegridEdges_t** self, size_t nVertsPerCell, size_t ncells, const double verts[]) {
    mnt_grid_new(&(*self)->dstGridObj);
    mnt_grid_setPointsPtr(&(*self)->dstGridObj, (int) nVertsPerCell, (vtkIdType) ncells, verts);
    mnt_grid_get(&(*self)->dstGridObj, &(*self)->dstGrid);
    return 0;
}

extern "C"
int mnt_regridedges_build(RegridEdges_t** self, int numCellsPerBucket) {

    // checks
    if (!(*self)->srcGrid) {
        std::cerr << "mnt_regridedges_build: ERROR must set source grid!\n";
        return 1;
    }
    if (!(*self)->dstGrid) {
        std::cerr << "mnt_regridedges_build: ERROR must set destination grid!\n";
        return 2;
    }

    // build the cell edge iterator
    QuadEdgeIter edgeIt;

    // build the locator
    (*self)->srcLoc->SetDataSet((*self)->srcGrid);
    (*self)->srcLoc->SetNumberOfCellsPerBucket(numCellsPerBucket);
    (*self)->srcLoc->BuildLocator(); 

    // compute the weights
    vtkIdList* dstPtIds = vtkIdList::New();
    vtkIdList* srcCellIds = vtkIdList::New();
    double dstEdgePt0[] = {0., 0., 0.};
    double dstEdgePt1[] = {0., 0., 0.};
    Vector<double> pcoords0(3);
    Vector<double> pcoords1(3);

    vtkPoints* dstPoints = (*self)->dstGrid->GetPoints();

    (*self)->numSrcCells = (*self)->srcGrid->GetNumberOfCells();
    (*self)->numDstCells = (*self)->dstGrid->GetNumberOfCells();

    // reserve some space for the weights and their cell/edge id arrays
    size_t n = (*self)->numDstCells * (*self)->numEdgesPerCell * 20;
    (*self)->weights.reserve(n);
    (*self)->weightSrcEdgeIds.reserve(n);
    (*self)->weightDstEdgeIds.reserve(n);
    (*self)->weightSrcCellIds.reserve(n);
    (*self)->weightDstCellIds.reserve(n);

    // iterate over the dst grid cells
    for (vtkIdType dstCellId = 0; dstCellId < (*self)->numDstCells; ++dstCellId) {

        // get this cell vertex Ids
        (*self)->dstGrid->GetCellPoints(dstCellId, dstPtIds);

        vtkCell* dstCell = (*self)->dstGrid->GetCell(dstCellId);
        int numEdges = dstCell->GetNumberOfEdges();

        // iterate over the four edges of each dst cell
        for (int dstEdgeIndex = 0; dstEdgeIndex < edgeIt.getNumberOfEdges(); ++dstEdgeIndex) {

            int id0, id1;
            edgeIt.getCellPointIds(dstEdgeIndex, &id0, &id1);
              
            dstPoints->GetPoint(dstCell->GetPointId(id0), dstEdgePt0);
            dstPoints->GetPoint(dstCell->GetPointId(id1), dstEdgePt1);

            // break the edge into sub-edges
            PolysegmentIter polySegIter = PolysegmentIter((*self)->srcGrid, 
                                                          (*self)->srcLoc,
                                                          dstEdgePt0, dstEdgePt1);

            // number of sub-segments
            size_t numSegs = polySegIter.getNumberOfSegments();

            // iterate over the sub-segments. Each sub-segment gets a src cell Id,
            // start/end cell param coords, the coefficient...
            polySegIter.reset();
            for (size_t iseg = 0; iseg < numSegs; ++iseg) {

                const vtkIdType srcCellId = polySegIter.getCellId();
                const Vector<double>& xia = polySegIter.getBegCellParamCoord();
                const Vector<double>& xib = polySegIter.getEndCellParamCoord();
                const double coeff = polySegIter.getCoefficient();

                Vector<double> dxi = xib - xia;
                Vector<double> xiMid = 0.5*(xia + xib);

                // create pair (dstCellId, srcCellId)
                std::pair<vtkIdType, vtkIdType> k = std::pair<vtkIdType, vtkIdType>(dstCellId, 
                                                                                    srcCellId);
                vtkCell* srcCell = (*self)->srcGrid->GetCell(srcCellId);
                double* srcCellParamCoords = srcCell->GetParametricCoords();

                for (int srcEdgeIndex = 0; srcEdgeIndex < edgeIt.getNumberOfEdges(); ++srcEdgeIndex) {

                    int i0, i1;
                    edgeIt.getCellPointIds(srcEdgeIndex, &i0, &i1);

                    // compute the interpolation weight, a product for every dimension
                    double weight = 1.0;
                    for (size_t d = 0; d < 2; ++d) { // 2d 

                        double xiM = xiMid[d];

                        // mid point of edge in parameter space
                        double x = 0.5*(srcCellParamCoords[i0*3 + d] + srcCellParamCoords[i1*3 + d]);

                        // use Lagrange interpolation to evaluate the basis function integral for
                        // any for the 3 possible x values in {0, 0.5, 1}. This formula will make 
                        // it easier to extend the code to 3d
                        double xm00 = x;
                        double xm05 = x - 0.5;
                        double xm10 = x - 1.0;
                        double lag00 = + 2 * xm05 * xm10;
                        double lag05 = - 4 * xm00 * xm10;
                        double lag10 = + 2 * xm00 * xm05;

                        // taking the abs value because the correct the sign for edges that 
                        // run from top to bottom or right to left.
                        weight *= (1.0 - xiM)*lag00 + dxi[d]*lag05 + xiM*lag10;
                    }

                    // coeff accounts for the duplicity in case where segments are shared between cells
                    weight *= coeff;

                    if (std::abs(weight) > 1.e-15) {
                        // only store the weights if they non-zero
                        (*self)->weights.push_back(weight);
                        (*self)->weightSrcCellIds.push_back(srcCellId);
                        (*self)->weightSrcEdgeIds.push_back(srcEdgeIndex);
                        (*self)->weightDstCellIds.push_back(dstCellId);
                        (*self)->weightDstEdgeIds.push_back(dstEdgeIndex);
                    }
                }

                // next segment
                polySegIter.next();

            }

            double totalT = polySegIter.getIntegratedParamCoord();
            if (std::abs(totalT - 1.0) > 1.e-10) {
            	printf("Warning: total t of segment: %lf != 1 (diff=%lg) dst cell %lld points (%18.16lf, %18.16lf), (%18.16lf, %18.16lf)\n",
            		   totalT, totalT - 1.0, dstCellId, dstEdgePt0[0], dstEdgePt0[1], dstEdgePt1[0], dstEdgePt1[1]);
            }

        }
    }

    // clean up
    srcCellIds->Delete();
    dstPtIds->Delete();

    return 0;
}

extern "C"
int mnt_regridedges_getNumSrcCells(RegridEdges_t** self, int* n) {
    *n = (*self)->numSrcCells;
    return 0;
}

extern "C"
int mnt_regridedges_getNumDstCells(RegridEdges_t** self, int* n) {
    *n = (*self)->numDstCells;
    return 0;
}

extern "C"
int mnt_regridedges_getNumEdgesPerCell(RegridEdges_t** self, int* n) {
    *n = (*self)->numEdgesPerCell;
    return 0;
}


extern "C"
int mnt_regridedges_applyWeights(RegridEdges_t** self, const double src_data[], double dst_data[]) {

    // initialize the data to zero
    size_t n = (*self)->numDstCells * (*self)->numEdgesPerCell;
    for (size_t i = 0; i < n; ++i) {
        dst_data[i] = 0.0;
    }

    // add the contributions from each cell overlaps
    for (size_t i = 0; i < (*self)->weights.size(); ++i) {

        vtkIdType dstCellId = (*self)->weightDstCellIds[i];
        vtkIdType srcCellId = (*self)->weightSrcCellIds[i];
        int dstEdgeIndex = (*self)->weightDstEdgeIds[i];
        int srcEdgeIndex = (*self)->weightSrcEdgeIds[i];

        // index into the flat array
        size_t dstK = dstEdgeIndex + (*self)->numEdgesPerCell * dstCellId;
        size_t srcK = srcEdgeIndex + (*self)->numEdgesPerCell * srcCellId;

        dst_data[dstK] += (*self)->weights[i] * src_data[srcK];
    }

    return 0;
}

extern "C"
int mnt_regridedges_load(RegridEdges_t** self, 
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

    // should check that numEdgesPerCell and (*self)->numEdgesPerCell match

    int dstCellIdsId, srcCellIdsId, dstEdgeIdsId, srcEdgeIdsId, weightsId;

    ier = nc_inq_varid(ncid, "dst_cell_ids", &dstCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }
    ier = nc_inq_varid(ncid, "src_cell_ids", &srcCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }
    ier = nc_inq_varid(ncid, "dst_edge_ids", &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 5;
    }
    ier = nc_inq_varid(ncid, "src_edge_ids", &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 6;
    }
    ier = nc_inq_varid(ncid, "weights", &weightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"weights\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 7;
    }

    (*self)->weights.resize(numWeights);
    (*self)->weightDstCellIds.resize(numWeights);
    (*self)->weightSrcCellIds.resize(numWeights);
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
    ier = nc_get_var_longlong(ncid, dstCellIdsId, &((*self)->weightDstCellIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not read var \"dst_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_get_var_longlong(ncid, srcCellIdsId, &((*self)->weightSrcCellIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_get_var_int(ncid, dstEdgeIdsId, &((*self)->weightDstEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 11;
    }
    ier = nc_get_var_int(ncid, srcEdgeIdsId, &((*self)->weightSrcEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 12;
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
int mnt_regridedges_dump(RegridEdges_t** self, 
		         const char* fort_filename, int n) {

    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);

    size_t numWeights = (*self)->weights.size();

    int ncid, ier;
    ier = nc_create(filename.c_str(), NC_CLOBBER, &ncid);
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

    int dstCellIdsId;
    ier = nc_def_var(ncid, "dst_cell_ids", NC_INT, 1, numWeightsAxis, &dstCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_cell_ids\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    int srcCellIdsId;
    ier = nc_def_var(ncid, "src_cell_ids", NC_INT, 1, numWeightsAxis, &srcCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    int dstEdgeIdsId;
    ier = nc_def_var(ncid, "dst_edge_ids", NC_INT, 1, numWeightsAxis, &dstEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 5;
    }

    int srcEdgeIdsId;
    ier = nc_def_var(ncid, "src_edge_ids", NC_INT, 1, numWeightsAxis, &srcEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 6;
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

    // write
    ier = nc_put_var_longlong(ncid, dstCellIdsId, &((*self)->weightDstCellIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"dst_cell_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }
    ier = nc_put_var_longlong(ncid, srcCellIdsId, &((*self)->weightSrcCellIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"src_cell_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_put_var_int(ncid, dstEdgeIdsId, &((*self)->weightDstEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"dst_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_put_var_int(ncid, srcEdgeIdsId, &((*self)->weightSrcEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"src_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 11;
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
int mnt_regridedges_print(RegridEdges_t** self) {
    size_t numWeights = (*self)->weights.size();
    QuadEdgeIter edgeIt;
    std::cout << "edge to vertex connectivity:\n";
    for (int edgeId = 0; edgeId < edgeIt.getNumberOfEdges(); ++edgeId) {
        int i0, i1;
        edgeIt.getCellPointIds(edgeId, &i0, &i1);
        std::cout << "edge " << edgeId << ": " << i0 << "->" << i1 << '\n';
    }
    std::cout << "Number of weights: " << numWeights << '\n';
    printf("                 dst_cell  dst_edge         src_cell  src_edge           weight\n");
    for (size_t i = 0; i < numWeights; ++i) {
    printf("%10ld       %8ld         %1d         %8ld         %1d   %15.5lg\n", i, 
               (*self)->weightDstCellIds[i], (*self)->weightDstEdgeIds[i], 
               (*self)->weightSrcCellIds[i], (*self)->weightSrcEdgeIds[i],
               (*self)->weights[i]);
    }
    return 0;
}

