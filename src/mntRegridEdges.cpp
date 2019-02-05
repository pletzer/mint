#include <mntRegridEdges.h>
#include <mntPolysegmentIter.h>
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
    (*self)->weights.clear();
    (*self)->numSrcCells = 0;
    (*self)->numDstCells = 0;
    (*self)->numPointsPerCell = 4; // 2d
    (*self)->numEdgesPerCell = 4; // 2d
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

    // iterate over the dst grid cells
    for (vtkIdType dstCellId = 0; dstCellId < (*self)->numDstCells; ++dstCellId) {

        // get this cell vertex Ids
        (*self)->dstGrid->GetCellPoints(dstCellId, dstPtIds);

        vtkCell* dstCell = (*self)->dstGrid->GetCell(dstCellId);
        int numEdges = dstCell->GetNumberOfEdges();

        // iterate over the four edges of each dst cell
        for (size_t i0 = 0; i0 < numEdges; ++i0) {

            vtkCell* dstEdge = dstCell->GetEdge(i0);
            vtkIdType id0 = dstEdge->GetPointId(0);
            vtkIdType id1 = dstEdge->GetPointId(1);
              
            dstPoints->GetPoint(id0, dstEdgePt0);
            dstPoints->GetPoint(id1, dstEdgePt1);

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

                // use vtkHexahedron in 3d!
                vtkQuad* srcQuad = dynamic_cast<vtkQuad*>(srcCell);

                // compute the weight contributions from each src edge, one value 
                // per edge
                std::vector<double> ws(numEdges, 1.0);
                        
                // iterate over the edges of the src cell
                for (size_t j0 = 0; j0 < numEdges; ++j0) {

                    // get the indices of the vertices for the edge
                    int* i01 = srcQuad->GetEdgeArray(j0);

                    // compute the interpolation weight, a product for every dimension
                    for (size_t d = 0; d < 3; ++d) {

                        double xiM = xiMid[d];
                        int i0 = i01[0];
                        int i1 = i01[1];

                        // mid point of edge in parameter space
                        double x = 0.5*(srcCellParamCoords[i0*3 + d] + srcCellParamCoords[i1*3 + d]);

                        // use Lagrange interpolation to evaluate the basis function integral for
                        // any for the 3 possible x values in {0, 0.5, 1}
                        double xm00 = x;
                        double xm05 = x - 0.5;
                        double xm10 = x - 1.0;
                        ws[j0] *= (1.0 - xiM) * xm05*xm10 \
                                -      dxi[d] * xm00*xm10 \
                                + (0.0 + xiM) * xm00*xm05;

                    }

                    // coeff accounts for the duplicity, some segments are shared between cells
                    ws[j0] *= 2*coeff; // WHY a coefficient 2? Comes from the Lagrange representation?
                }

                if ((*self)->weights.find(k) == (*self)->weights.end()) {
                    // no current entry for dstCellId and srcCellIds
                    // initialize the weights to a zero matrix
                    std::vector<double> zeros( (*self)->numEdgesPerCell, 0.0 );
                    std::pair< std::pair<vtkIdType, vtkIdType>, std::vector<double> > kv 
                    // create an entry (dstCellId, srcCellId) -> zeros(numEdges)
                      = std::pair< std::pair<vtkIdType, vtkIdType>, std::vector<double> >(k, zeros);
                    (*self)->weights.insert(kv);
                }

                // search for entry
                std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> >::iterator 
                     it = (*self)->weights.find(k);

                // add the weights
                for (size_t j = 0; j < it->second.size(); ++j) {
                    it->second[j] += ws[j];
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
    size_t numEdges = (*self)->numDstCells * (*self)->numEdgesPerCell;
    for (size_t i = 0; i < numEdges; ++i) {
        dst_data[i] = 0.0;
    }

    for (std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> >::const_iterator 
         it = (*self)->weights.begin(); it != (*self)->weights.end(); ++it) {

        vtkIdType dstCellId = it->first.first;
        vtkIdType srcCellId = it->first.second;
        const std::vector<double>& weights = it->second;

        size_t kd = dstCellId * (*self)->numEdgesPerCell;
        size_t ks = srcCellId * (*self)->numEdgesPerCell;

        for (size_t ie = 0; ie < (*self)->numEdgesPerCell; ++ie) {
            dst_data[kd + ie] += weights[ie] * src_data[ks + ie];
        }
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

    // get the sizes
    size_t numWeights, numEdgesPerCell;
    int numWeightsId;
    int numEdgesId;
    ier = nc_inq_dimid(ncid, "num_weights", &numWeightsId);
    ier = nc_inq_dimlen(ncid, numWeightsId, &numWeights);
    ier = nc_inq_dimid(ncid, "num_edges_per_cell", &numEdgesId);
    ier = nc_inq_dimlen(ncid, numEdgesId, &numEdgesPerCell);

    // should check that numEdgesPerCell and (*self)->numEdgesPerCell match

    int dstCellIdsId, srcCellIdsId, weightsId;

    ier = nc_inq_varid(ncid, "dst_cell_ids", &dstCellIdsId);
    ier = nc_inq_varid(ncid, "src_cell_ids", &srcCellIdsId);
    ier = nc_inq_varid(ncid, "weights", &weightsId);

    // read
    std::vector<long> dstCellIds(numWeights);
    std::vector<long> srcCellIds(numWeights);
    std::vector<double> weights(numWeights * numEdgesPerCell);
    ier = nc_get_var_long(ncid, dstCellIdsId, &dstCellIds[0]);
    ier = nc_get_var_long(ncid, srcCellIdsId, &srcCellIds[0]);
    ier = nc_get_var_double(ncid, weightsId, &weights[0]);

    ier = nc_close(ncid);    

    // store in map
    for (int i = 0; i < numWeights; ++i) {

        // create (dstCellId, srcCellId) entry
        std::pair< vtkIdType, vtkIdType > k(dstCellIds[i], srcCellIds[i]);

        // create vector of weights for this cell
        std::vector<double> v(&weights[numEdgesPerCell*i], 
                              &weights[numEdgesPerCell*i + numEdgesPerCell]);

        // create pair of entries
        std::pair< std::pair<vtkIdType, vtkIdType>, std::vector<double> > kv(k, v);

        // insert
        (*self)->weights.insert(kv);
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

    // create dimensions
    int numWeightsId;
    int numEdgesId;
    ier = nc_def_dim(ncid, "num_weights", (int) numWeights, &numWeightsId);
    ier = nc_def_dim(ncid, "num_edges_per_cell", (*self)->numEdgesPerCell, &numEdgesId);


    // create variables
    int numWeightsAxis[] = {numWeightsId};
    int numWeightsNumEdgesAxes[] = {numWeightsId, numEdgesId};

    int dstCellIdsId;
    ier = nc_def_var(ncid, "dst_cell_ids", NC_LONG, 1, numWeightsAxis, &dstCellIdsId);

    int srcCellIdsId;
    ier = nc_def_var(ncid, "src_cell_ids", NC_LONG, 1, numWeightsAxis, &srcCellIdsId);

    int weightsId;
    ier = nc_def_var(ncid, "weights", NC_DOUBLE, 2, numWeightsNumEdgesAxes, &weightsId);

    // close define mode
    ier = nc_enddef(ncid);

    // load into arrays
    std::vector<long> dstCellIds(numWeights);
    std::vector<long> srcCellIds(numWeights);
    std::vector<double> weights(numWeights * (*self)->numEdgesPerCell);

    // dst cell, src cell counter
    size_t i = 0;
    for (std::map< std::pair<vtkIdType, vtkIdType>, std::vector<double> >::const_iterator 
        it = (*self)->weights.begin(); it != (*self)->weights.end(); ++it) {

        // get the dst/src cell indices
        dstCellIds[i] = it->first.first;
        srcCellIds[i] = it->first.second;

        // iterate over the edges
        for (int j = 0; j < (*self)->numEdgesPerCell; ++j) {
            weights[(*self)->numEdgesPerCell*i + j] = it->second[j];
        }
        i++;
    }

    // write
    ier = nc_put_var_long(ncid, dstCellIdsId, &dstCellIds[0]);
    ier = nc_put_var_long(ncid, srcCellIdsId, &srcCellIds[0]);
    ier = nc_put_var_double(ncid, weightsId, &weights[0]);

    ier = nc_close(ncid);

    return 0;
}


