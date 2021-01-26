#include <mntRegridEdges.h>
#include <mntPolysegmentIter.h>
#include <mntNcFieldRead.h>
#include <mntNcFieldWrite.h>
#include <mntWeights.h>

#include <netcdf.h>

#include <iostream>
#include <cstdio>
#include <cstring>

#include <vtkIdList.h>
#include <vtkHexahedron.h> // for 3d grids
#include <vtkQuad.h>       // for 2d grids
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPoints.h>


extern "C"
int mnt_regridedges_new(RegridEdges_t** self) {
    
    *self = new RegridEdges_t();
    (*self)->srcGrid = NULL;
    (*self)->dstGrid = NULL;
    (*self)->srcLoc = vmtCellLocator::New();
    (*self)->numPointsPerCell = 4; // 2d
    (*self)->numEdgesPerCell = 4;  // 2d

    mnt_grid_new(&((*self)->srcGridObj));
    mnt_grid_new(&((*self)->dstGridObj));

    (*self)->ndims = 0;
    // multiarray iterator
    (*self)->mai = NULL;

    (*self)->srcReader = NULL;
    (*self)->srcNcid = -1;
    (*self)->srcVarid = -1;

    (*self)->dstWriter = NULL;

    return 0;
}

extern "C"
int mnt_regridedges_del(RegridEdges_t** self) {

    int ier = 0;

    // destroy the cell locator
    (*self)->srcLoc->Delete();
   
    // destroy the source and destination grids
    mnt_grid_del(&((*self)->srcGridObj));
    mnt_grid_del(&((*self)->dstGridObj));

    if ((*self)->srcNcid >= 0) {
        ier = nc_close((*self)->srcNcid);
    }

    if ((*self)->mai) {
        ier = mnt_multiarrayiter_del(&(*self)->mai);
    }

    if ((*self)->srcReader) {
        ier = mnt_ncfieldread_del(&(*self)->srcReader);
    }

    if ((*self)->dstWriter) {
        ier = mnt_ncfieldwrite_del(&(*self)->dstWriter);
    }

    delete *self;

    return ier;
}

extern "C"
int mnt_regridedges_setSrcGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole) {

    return mnt_grid_setFlags( &(*self)->srcGridObj, fixLonAcrossDateline, averageLonAtPole );

}

extern "C"
int mnt_regridedges_setDstGridFlags(RegridEdges_t** self, int fixLonAcrossDateline, int averageLonAtPole) {

    return mnt_grid_setFlags( &(*self)->dstGridObj, fixLonAcrossDateline, averageLonAtPole );

}

extern "C"
int mnt_regridedges_dumpSrcGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength) {

    return mnt_grid_dump( &(*self)->srcGridObj, std::string(fort_filename).c_str() );

}

extern "C"
int mnt_regridedges_dumpDstGridVtk(RegridEdges_t** self,
                                   const char* fort_filename, int nFilenameLength) {

    return mnt_grid_dump( &(*self)->dstGridObj, std::string(fort_filename).c_str() );

}


extern "C"
int mnt_regridedges_initSliceIter(RegridEdges_t** self,
                                  const char* src_fort_filename, int src_nFilenameLength,
                                  const char* dst_fort_filename, int dst_nFilenameLength,
                                  int append,
                                  const char* field_name, int nFieldNameLength, 
                                  size_t* numSlices) {

    int ier = 0;

    std::string srcFileAndMeshName = std::string(src_fort_filename, src_nFilenameLength);
    std::string dstFileAndMeshName = std::string(dst_fort_filename, dst_nFilenameLength);
    std::string fieldname = std::string(field_name, nFieldNameLength);

    // filter out the mesh name, if present (not used here)
    size_t columnL = srcFileAndMeshName.find(':');
    std::string srcFilename = srcFileAndMeshName.substr(0, columnL);
    columnL = dstFileAndMeshName.find(':');
    std::string dstFilename = dstFileAndMeshName.substr(0, columnL);

    // open the source file
    ier = nc_open(srcFilename.c_str(), NC_NOWRITE, &(*self)->srcNcid);
    if (ier != 0) {
        std::cerr << "ERROR: could not open " << srcFilename << '\n';
        return 1;
    }

    // get tht variable id
    ier = nc_inq_varid((*self)->srcNcid, fieldname.c_str(), &(*self)->srcVarid);
    if (ier != 0) {
        std::cerr << "ERROR: could not find variable " << fieldname << " in file " << srcFilename << '\n';
        return 1;
    }

    // intantiate the reader
    ier = mnt_ncfieldread_new(&(*self)->srcReader, (*self)->srcNcid, (*self)->srcVarid);

    // get the number of dimensions
    ier = mnt_ncfieldread_getNumDims(&(*self)->srcReader, &(*self)->ndims);
    if (ier != 0) {
        std::cerr << "ERROR: getting the number of dims of " << fieldname << " from file " << srcFilename << '\n';
        return 2;
    }

    // allocate
    (*self)->startIndices.resize((*self)->ndims, 0); // initialize to zero
    (*self)->dimNames.resize((*self)->ndims);
    // slice has dimension one except for the edge axis
    (*self)->srcCounts.resize((*self)->ndims, 1);
    (*self)->dstCounts.resize((*self)->ndims, 1); 
    (*self)->srcDims.resize((*self)->ndims, 0);
    (*self)->dstDims.resize((*self)->ndims, 0);

    //
    // assume that the src and dst data have the same axes/dimensions except for the last (number of edges)
    //
    for (int i = 0; i < (*self)->ndims - 1; ++i) {
        // get the source field dimensions along each axis
        ier = mnt_ncfieldread_getDim(&(*self)->srcReader, i, &(*self)->srcDims[i]);
        if (ier != 0) {
            std::cerr << "ERROR: getting the dimension " << i << " from source file\n";
        }

        
        (*self)->dimNames[i].resize(256);


        ier = mnt_ncfieldread_getDimName(&(*self)->srcReader, i, 
                                         &(*self)->dimNames[i][0], (int) (*self)->dimNames[i].size());
        if (ier != 0) {
            std::cerr << "ERROR: getting the dimension name " << i << " from source file\n";
        }

        // all except last dimensions are the same 
        (*self)->dstDims[i] = (*self)->srcDims[i];
    }

    // last dimension is edge axis
    int i = (*self)->ndims - 1;
    size_t numSrcEdges;
    size_t numDstEdges;
    ier = mnt_grid_getNumberOfEdges(&(*self)->srcGridObj, &numSrcEdges);
    ier = mnt_grid_getNumberOfEdges(&(*self)->dstGridObj, &numDstEdges);
    (*self)->srcCounts[i] = numSrcEdges;
    (*self)->dstCounts[i] = numDstEdges;
    (*self)->srcDims[i] = numSrcEdges;
    (*self)->dstDims[i] = numDstEdges;
    (*self)->dimNames[i].resize(128);
    ier = mnt_ncfieldread_getDimName(&(*self)->srcReader, i, 
                                     &(*self)->dimNames[i][0], (int) (*self)->dimNames[i].size());

    // initialize the writer
    ier = mnt_ncfieldwrite_new(&(*self)->dstWriter, dstFilename.c_str(), (int) dstFilename.size(), 
                                fieldname.c_str(), (int) fieldname.size(), append);
    if (ier != 0) {
        std::cerr << "ERROR: occurred when creating/opening file " << dstFilename << " with field " 
                  << fieldname << " in append mode " << append << '\n';
        return 1;
    }

    ier = mnt_ncfieldwrite_setNumDims(&(*self)->dstWriter, (*self)->ndims);
    if (ier != 0) {
        std::cerr << "ERROR: cannot set the number of dimensions for field " << fieldname << " in file " << dstFilename << '\n';
        return 2;
    }

    // add num_edges axis
    for (int i = 0; i < (*self)->ndims; ++i) {
        std::string axname = (*self)->dimNames[i];
        ier = mnt_ncfieldwrite_setDim(&(*self)->dstWriter, i, axname.c_str(), (int) axname.size(), (*self)->dstDims[i]);
        if (ier != 0) {
            std::cerr << "ERROR: setting dimension " << i << " (" << axname << ") to " << (*self)->dstDims[i]
                  << " for field " << fieldname << " in file " << dstFilename << '\n';
            return 3;
        }
    }

    // create iterator, assume the last dimension is the number of edges. Note ndims - 1
    ier = mnt_multiarrayiter_new(&(*self)->mai, (*self)->ndims - 1, &(*self)->srcDims[0]);
    if (ier != 0) {
        std::cerr << "ERROR: creating multiarray iterator\n";
        return 4;
    }

    ier = mnt_multiarrayiter_getNumIters(&(*self)->mai, numSlices);
    if (ier != 0) {
        std::cerr << "ERROR: getting the number of iterations from the multiarray iterator\n";
        return 4;
    }

    return ier;
}


extern "C"
int mnt_regridedges_loadSrcSlice(RegridEdges_t** self,
                                 double data[]) {

    if (!(*self)->srcReader) {
        std::cerr << "ERROR: must call mnt_regridedges_initSliceIter prior to mnt_regridedges_loadSrcSlice\n";
        return 5;        
    }

    // get the current slice indices
    int ier = mnt_multiarrayiter_getIndices(&(*self)->mai, &(*self)->startIndices[0]);


    ier = mnt_ncfieldread_dataSlice(&(*self)->srcReader, 
                                    &(*self)->startIndices[0], 
                                    &(*self)->srcCounts[0], data);
    if (ier != 0) {
        std::cerr << "ERROR: occurred when loading slice of src data\n";
        return 4;
    }

    return 0;
}

extern "C"
int mnt_regridedges_dumpDstSlice(RegridEdges_t** self,
                                 double data[]) {

    if (!(*self)->dstWriter) {
        std::cerr << "ERROR: must call mnt_regridedges_initSliceIter prior to mnt_regridedges_dumpDstSlice\n";
        return 5;        
    }

    // get the current slice indices
    int ier = mnt_multiarrayiter_getIndices(&(*self)->mai, &(*self)->startIndices[0]);


    ier = mnt_ncfieldwrite_dataSlice(&(*self)->dstWriter, 
                                     &(*self)->startIndices[0], 
                                     &(*self)->dstCounts[0], data);
    if (ier != 0) {
        std::cerr << "ERROR: occurred when dumping slice of dst data\n";
        return 4;
    }

    return 0;
}


extern "C"
int mnt_regridedges_nextSlice(RegridEdges_t** self) {

    // increment the iterator
    int ier = mnt_multiarrayiter_next(&(*self)->mai);
    
    return ier;
}


extern "C"
int mnt_regridedges_loadEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, double data[]) {

    int ier = 0;

    std::string fileAndMeshName = std::string(fort_filename, nFilenameLength);

    // filter out the mesh name, if present (not used here)
    size_t columnL = fileAndMeshName.find(':');
    std::string filename = fileAndMeshName.substr(0, columnL);

    std::string fieldname = std::string(field_name, nFieldNameLength);

    int ncid;
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != 0) {
        std::cerr << "ERROR: could not open " << filename << '\n';
        return 1;
    }

    int varid;
    ier = nc_inq_varid(ncid, fieldname.c_str(), &varid);
    if (ier != 0) {
        std::cerr << "ERROR: could not find variable " << fieldname << " in file " << filename << '\n';
        nc_close(ncid);
        return 1;
    }

    NcFieldRead_t* rd = NULL;
    ier = mnt_ncfieldread_new(&rd, ncid, varid);

    // get the number of dimensions
    int ndims;
    ier = mnt_ncfieldread_getNumDims(&rd, &ndims);
    if (ier != 0) {
        std::cerr << "ERROR: getting the number of dims of " << fieldname << " from file " << filename << '\n';
        ier = mnt_ncfieldread_del(&rd);
        nc_close(ncid);
        return 2;
    }

    if (ndims != 1) {
        std::cerr << "ERROR: number of dimensions must be 1, got " << ndims << '\n';
        ier = mnt_ncfieldread_del(&rd);
        nc_close(ncid);
        return 3;        
    }

    ier = mnt_ncfieldread_data(&rd, data);
    if (ier != 0) {
        std::cerr << "ERROR: reading field " << fieldname << " from file " << filename << '\n';
        ier = mnt_ncfieldread_del(&rd);
        nc_close(ncid);
        return 4;
    }

    ier = mnt_ncfieldread_del(&rd);
    nc_close(ncid);

    return 0;
}


extern "C"
int mnt_regridedges_dumpEdgeField(RegridEdges_t** self,
                                  const char* fort_filename, int nFilenameLength,
                                  const char* field_name, int nFieldNameLength,
                                  size_t ndata, const double data[]) {
    
    std::string fileAndMeshName = std::string(fort_filename, nFilenameLength);
    std::string fieldname = std::string(field_name, nFieldNameLength);

    size_t columnL = fileAndMeshName.find(':');

    // get the file name
    std::string filename = fileAndMeshName.substr(0, columnL);
    // get the mesh name
    std::string meshname = fileAndMeshName.substr(columnL + 1);

    int ier;
    NcFieldWrite_t* wr = NULL;

    int n1 = filename.size();
    int n2 = fieldname.size();
    const int append = 0; // new file
    ier = mnt_ncfieldwrite_new(&wr, filename.c_str(), n1, fieldname.c_str(), n2, append);
    if (ier != 0) {
        std::cerr << "ERROR: create file " << filename << " with field " 
                  << fieldname << " in append mode " << append << '\n';
        return 1;
    }

    ier = mnt_ncfieldwrite_setNumDims(&wr, 1); // 1D array only in this implementation
    if (ier != 0) {
        std::cerr << "ERROR: cannot set the number of dimensions for field " << fieldname << " in file " << filename << '\n';
        ier = mnt_ncfieldwrite_del(&wr);
        return 2;
    }

    // add num_edges axis
    std::string axname = "num_edges";
    int n3 = axname.size();
    ier = mnt_ncfieldwrite_setDim(&wr, 0, axname.c_str(), n3, ndata);
    if (ier != 0) {
        std::cerr << "ERROR: setting dimension 0 (" << axname << ") to " << ndata 
                  << " for field " << fieldname << " in file " << filename << '\n';
        ier = mnt_ncfieldwrite_del(&wr);
        return 3;
    }

    // write the data to disk
    ier = mnt_ncfieldwrite_data(&wr, data);
    if (ier != 0) {
        std::cerr << "ERROR: writing data for field " << fieldname << " in file " << filename << '\n';
        ier = mnt_ncfieldwrite_del(&wr);
        return 5;
    }

    // clean up
    ier = mnt_ncfieldwrite_del(&wr);

    return 0;
}


extern "C"
int mnt_regridedges_loadSrcGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    int ier = mnt_grid_loadFrom2DUgrid(&((*self)->srcGridObj), filename.c_str());
    (*self)->srcGrid = (*self)->srcGridObj->grid;
    return ier;
}

extern "C"
int mnt_regridedges_loadDstGrid(RegridEdges_t** self, 
                                const char* fort_filename, int n) {
    // Fortran strings don't come with null-termination character. Copy string 
    // into a new one and add '\0'
    std::string filename = std::string(fort_filename, n);
    int ier = mnt_grid_loadFrom2DUgrid(&((*self)->dstGridObj), filename.c_str());
    (*self)->dstGrid = (*self)->dstGridObj->grid;
    return ier;
}

extern "C"
int mnt_regridedges_build(RegridEdges_t** self, int numCellsPerBucket, double periodX, int debug) {

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
    (*self)->srcLoc->setPeriodicityLengthX(periodX);

    // compute the weights
    vtkIdList* dstPtIds = vtkIdList::New();
    vtkIdList* srcCellIds = vtkIdList::New();
    double dstEdgePt0[] = {0., 0., 0.};
    double dstEdgePt1[] = {0., 0., 0.};
    vtkPoints* dstPoints = (*self)->dstGrid->GetPoints();

    size_t numDstCells = (*self)->dstGrid->GetNumberOfCells();

    // reserve some space for the weights and their cell/edge id arrays
    size_t n = numDstCells * (*self)->numEdgesPerCell * 20;
    (*self)->weights.reserve(n);
    (*self)->weightSrcFaceEdgeIds.reserve(n);
    (*self)->weightDstFaceEdgeIds.reserve(n);
    (*self)->weightSrcCellIds.reserve(n);
    (*self)->weightDstCellIds.reserve(n);

    vtkPoints* badSegmentsPoints = NULL;
    vtkUnstructuredGrid* badSegmentsGrid = NULL;
    vtkIdList* badSegmentPtIds = NULL;
    vtkIdType badPtId = 0;
    if (debug == 3) {
        printf("   dstCellId dstEdgeIndex     dstEdgePt0     dstEdgePt1     srcCellId            xia          xib        ta     tb   tmax\n");
    }
    else if (debug == 2) {
        badSegmentsPoints = vtkPoints::New();
        badSegmentsGrid = vtkUnstructuredGrid::New();
        badSegmentsGrid->SetPoints(badSegmentsPoints);
        badSegmentPtIds = vtkIdList::New();
        badSegmentPtIds->SetNumberOfIds(2);
    }

    // iterate over the dst grid cells
    int numBadSegments = 0;
    for (size_t dstCellId = 0; dstCellId < numDstCells; ++dstCellId) {

        // get this cell vertex Ids
        (*self)->dstGrid->GetCellPoints(dstCellId, dstPtIds);

        vtkCell* dstCell = (*self)->dstGrid->GetCell(dstCellId);

        // iterate over the four edges of each dst cell
        for (int dstEdgeIndex = 0; dstEdgeIndex < (*self)->edgeConnectivity.getNumberOfEdges(); 
             ++dstEdgeIndex) {

            int id0, id1;
            (*self)->edgeConnectivity.getCellPointIds(dstEdgeIndex, &id0, &id1);
            
            // fill in the start/end points of this edge  
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
                const Vec3& xia = polySegIter.getBegCellParamCoord();
                const Vec3& xib = polySegIter.getEndCellParamCoord();
                const double coeff = polySegIter.getCoefficient();

                if (debug == 3) {
                    printf("%12zu %12d    %5.3lf,%5.3lf    %5.3lf,%5.3lf  %12lld    %5.3lf,%5.3lf  %5.3lf,%5.3lf   %5.4lf, %5.4lf   %10.7lf\n", 
                        dstCellId, dstEdgeIndex, 
                        dstEdgePt0[0], dstEdgePt0[1], 
                        dstEdgePt1[0], dstEdgePt1[1], 
                        srcCellId,
                        xia[0], xia[1], xib[0], xib[1], 
                        polySegIter.getBegLineParamCoord(), polySegIter.getEndLineParamCoord(),
                        polySegIter.getIntegratedParamCoord());
                }

                vtkCell* srcCell = (*self)->srcGrid->GetCell(srcCellId);
                double* srcCellParamCoords = srcCell->GetParametricCoords();

                for (int srcEdgeIndex = 0; srcEdgeIndex < (*self)->edgeConnectivity.getNumberOfEdges(); 
                       ++srcEdgeIndex) {

                    int is0, is1;
                    (*self)->edgeConnectivity.getCellPointIds(srcEdgeIndex, &is0, &is1);
                    
                    // compute the interpolation weight
                    double weight = computeWeight(&srcCellParamCoords[is0*3], 
                                                  &srcCellParamCoords[is1*3], 
                                                  xia, xib);

                    // coeff accounts for the duplicity in the case where segments are shared between cells
                    weight *= coeff;

                    if (std::abs(weight) > 1.e-15) {
                        // only store the weights if they non-zero
                        (*self)->weights.push_back(weight);
                        (*self)->weightSrcCellIds.push_back(srcCellId);
                        (*self)->weightSrcFaceEdgeIds.push_back(srcEdgeIndex);
                        (*self)->weightDstCellIds.push_back(dstCellId);
                        (*self)->weightDstFaceEdgeIds.push_back(dstEdgeIndex);
                    }
                }

                // next segment
                polySegIter.next();

            }

            if (debug > 0) {
                double totalT = polySegIter.getIntegratedParamCoord();
                if (std::abs(totalT - 1.0) > 1.e-10) {
                    printf("Warning: [%d] total t of segment: %lf != 1 (diff=%lg) dst cell %zu points (%18.16lf, %18.16lf), (%18.16lf, %18.16lf)\n",
                       numBadSegments, totalT, totalT - 1.0, dstCellId, dstEdgePt0[0], dstEdgePt0[1], dstEdgePt1[0], dstEdgePt1[1]);
                    numBadSegments++;

                    if (debug == 2) {
                        badSegmentsPoints->InsertNextPoint(dstEdgePt0);
                        badSegmentsPoints->InsertNextPoint(dstEdgePt1);
                        badSegmentPtIds->SetId(0, badPtId);
                        badSegmentPtIds->SetId(1, badPtId + 1);
                        badSegmentsGrid->InsertNextCell(VTK_LINE, badSegmentPtIds);
                        badPtId += 2;
                    }
                }
            }

        }
    }

    // clean up
    srcCellIds->Delete();
    dstPtIds->Delete();

    if (debug == 2 && badPtId > 0) {
        vtkUnstructuredGridWriter* wr = vtkUnstructuredGridWriter::New();
        std::string fname = "badSegments.vtk";
        std::cout << "Warning: saving segments that are not fully contained in the source grid in file " << fname << '\n';
        wr->SetFileName(fname.c_str());
        wr->SetInputData(badSegmentsGrid);
        wr->Update();
        wr->Delete();
        badSegmentPtIds->Delete();
        badSegmentsGrid->Delete();
        badSegmentsPoints->Delete();
    }

    return 0;
}

extern "C"
int mnt_regridedges_getNumSrcCells(RegridEdges_t** self, size_t* n) {
    *n = (*self)->srcGrid->GetNumberOfCells();
    return 0;
}

extern "C"
int mnt_regridedges_getNumDstCells(RegridEdges_t** self, size_t* n) {
    *n = (*self)->dstGrid->GetNumberOfCells();
    return 0;
}

extern "C"
int mnt_regridedges_getNumEdgesPerCell(RegridEdges_t** self, int* n) {
    *n = (*self)->numEdgesPerCell;
    return 0;
}

extern "C"
int mnt_regridedges_getNumSrcEdges(RegridEdges_t** self, size_t* nPtr) {
    if (!(*self)->srcGridObj) {
        std::cerr << "ERROR: source grid was not loaded\n";
        return 1;
    }
    int ier = mnt_grid_getNumberOfEdges(&((*self)->srcGridObj), nPtr);
    return ier;
}

extern "C"
int mnt_regridedges_getNumDstEdges(RegridEdges_t** self, size_t* nPtr) {
    if (!(*self)->dstGridObj) {
        std::cerr << "ERROR: destination grid was not loaded\n";
        return 1;
    }
    int ier = mnt_grid_getNumberOfEdges(&((*self)->dstGridObj), nPtr);
    return ier;
}

extern "C"
int mnt_regridedges_apply(RegridEdges_t** self, 
                          const double src_data[], double dst_data[]) {


    // make sure (*self)->srcGridObj.faceNodeConnectivity and the rest have been allocated
    if (!(*self)->srcGridObj ||
        (*self)->srcGridObj->faceNodeConnectivity.size() == 0 || 
        (*self)->srcGridObj->faceEdgeConnectivity.size() == 0 ||
        (*self)->srcGridObj->edgeNodeConnectivity.size() == 0) {
        std::cerr << "ERROR: looks like the src grid connectivity is not set.\n";
        std::cerr << "Typically this would occur if you did not read the grid\n";
        std::cerr << "from the netcdf Ugrid file.\n";
        return 1;
    }

    int ier;

    // number of unique edges on the destination grid
    size_t numDstEdges;
    ier = mnt_grid_getNumberOfEdges(&((*self)->dstGridObj), &numDstEdges);
    
    // initialize the data to zero
    for (size_t i = 0; i < numDstEdges; ++i) {
        dst_data[i] = 0.0;
    }

    // number of faces sharing the same edge. As a result of the multivaluedness of the longitudes,
    // each of the edges is treated independently (the longitude values may be different as seen from
    // the two faces). edgeMultiplicity tracks the number of adjacent faces. 
    std::vector<int> edgeMultiplicity(numDstEdges, 0);

    // add the contributions from each cell overlaps
    for (size_t i = 0; i < (*self)->weights.size(); ++i) {

        vtkIdType dstCellId = (*self)->weightDstCellIds[i];
        vtkIdType srcCellId = (*self)->weightSrcCellIds[i];
        int dstEdgeIndex = (*self)->weightDstFaceEdgeIds[i];
        int srcEdgeIndex = (*self)->weightSrcFaceEdgeIds[i];

        size_t srcEdgeId, dstEdgeId;
        int srcEdgeSign, dstEdgeSign;
        ier = mnt_grid_getEdgeId(&((*self)->srcGridObj), srcCellId, srcEdgeIndex, &srcEdgeId, &srcEdgeSign);
        ier = mnt_grid_getEdgeId(&((*self)->dstGridObj), dstCellId, dstEdgeIndex, &dstEdgeId, &dstEdgeSign);

        dst_data[dstEdgeId] += srcEdgeSign * dstEdgeSign * (*self)->weights[i] * src_data[srcEdgeId];

        // up to 2 faces can own this edge
        edgeMultiplicity[dstEdgeId] = std::min(2, edgeMultiplicity[dstEdgeId] + 1);
    }

    for (size_t i = 0; i < numDstEdges; ++i) {

        // there has been cases where edgeMultiplicity[i] is zero and so we need to guard 
        // against a division by zero. I would expect in this case dst_data[i] to be also 
        // zero but this would need to be checked! (edgeMultiplicity[i] is zero if the dst 
        // edge lies outside the src domain)
        dst_data[i] /= std::max(1, edgeMultiplicity[i]);
    }

    return 0;
}


extern "C"
int mnt_regridedges_loadWeights(RegridEdges_t** self, 
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

    int dstCellIdsId, srcCellIdsId, dstFaceEdgeIdsId, srcFaceEdgeIdsId, weightsId;

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
    ier = nc_inq_varid(ncid, "dst_face_edge_ids", &dstFaceEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_face_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 5;
    }
    ier = nc_inq_varid(ncid, "src_face_edge_ids", &srcFaceEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_face_edge_ids\"!\n";
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
    (*self)->weightDstFaceEdgeIds.resize(numWeights);
    (*self)->weightSrcFaceEdgeIds.resize(numWeights);

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
    ier = nc_get_var_int(ncid, dstFaceEdgeIdsId, &((*self)->weightDstFaceEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"dst_face_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 11;
    }
    ier = nc_get_var_int(ncid, srcFaceEdgeIdsId, &((*self)->weightSrcFaceEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not get ID for var \"src_face_edge_ids\"!\n";
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
int mnt_regridedges_dumpWeights(RegridEdges_t** self, 
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

    int numSpaceDimsId;
    ier = nc_def_dim(ncid, "num_space_dims", 3, &numSpaceDimsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_space_dims\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }    

    int numEdgesPerCellId;
    ier = nc_def_dim(ncid, "num_edges_per_cell", (*self)->numEdgesPerCell, &numEdgesPerCellId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_edges_per_cell\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }

    int numWeightsId;
    ier = nc_def_dim(ncid, "num_weights", (int) numWeights, &numWeightsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define dimension \"num_weights\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 2;
    }

    // create variables
    int paramCoordsAxis[] = {numEdgesPerCellId, numSpaceDimsId};
    int numWeightsAxis[] = {numWeightsId};

    int edgeParamCoordBegId, edgeParamCoordEndId;
    ier = nc_def_var(ncid, "edge_param_coord_beg", NC_DOUBLE, 2, paramCoordsAxis, &edgeParamCoordBegId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"edge_param_coord_beg\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }
    ier = nc_def_var(ncid, "edge_param_coord_end", NC_DOUBLE, 2, paramCoordsAxis, &edgeParamCoordEndId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"edge_param_coord_end\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    int dstCellIdsId;
    ier = nc_def_var(ncid, "dst_cell_ids", NC_INT64, 1, numWeightsAxis, &dstCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_cell_ids\"! ier = " << ier << "\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 3;
    }

    int srcCellIdsId;
    ier = nc_def_var(ncid, "src_cell_ids", NC_INT64, 1, numWeightsAxis, &srcCellIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_cell_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 4;
    }

    int dstFaceEdgeIdsId;
    ier = nc_def_var(ncid, "dst_face_edge_ids", NC_INT, 1, numWeightsAxis, &dstFaceEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"dst_face_edge_ids\"!\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 5;
    }

    int srcFaceEdgeIdsId;
    ier = nc_def_var(ncid, "src_face_edge_ids", NC_INT, 1, numWeightsAxis, &srcFaceEdgeIdsId);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not define variable \"src_face_edge_ids\"!\n";
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

    double edgeParamCoordBegs[(*self)->numEdgesPerCell * 3];
    double edgeParamCoordEnds[(*self)->numEdgesPerCell * 3];
    double* xiBeg;
    double* xiEnd;
    for (size_t e = 0; e < (*self)->numEdgesPerCell; ++e) {
        (*self)->edgeConnectivity.getParamCoords(e, &xiBeg, &xiEnd);
        for (size_t j = 0; j < 3; ++j) { // always 3d
            edgeParamCoordBegs[e*3 + j] = xiBeg[j];
            edgeParamCoordEnds[e*3 + j] = xiEnd[j];
        }
    }

    // write
    ier = nc_put_var_double(ncid, edgeParamCoordBegId, edgeParamCoordBegs);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"edge_param_coord_beg\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }

    ier = nc_put_var_double(ncid, edgeParamCoordEndId, edgeParamCoordEnds);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"edge_param_coord_end\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 9;
    }

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
    ier = nc_put_var_int(ncid, dstFaceEdgeIdsId, &((*self)->weightDstFaceEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"dst_face_edge_ids\"\n";
        std::cerr << nc_strerror (ier);
        nc_close(ncid);
        return 10;
    }
    ier = nc_put_var_int(ncid, srcFaceEdgeIdsId, &((*self)->weightSrcFaceEdgeIds)[0]);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: could not write variable \"src_face_edge_ids\"\n";
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
    std::cout << "edge to vertex connectivity:\n";
    for (int faceEdgeId = 0; 
         faceEdgeId < (*self)->edgeConnectivity.getNumberOfEdges(); 
         ++faceEdgeId) {
        int i0, i1;
        (*self)->edgeConnectivity.getCellPointIds(faceEdgeId, &i0, &i1);
        std::cout << "edge " << faceEdgeId << ": " << i0 << "->" << i1 << '\n';
    }
    std::cout << "Number of weights: " << numWeights << '\n';
    printf("                 dst_cell  dst_face_edge     src_cell  src_face_edge       weight\n");
    for (size_t i = 0; i < numWeights; ++i) {
    printf("%10ld       %8lld         %1d          %8lld         %1d   %15.5lg\n", 
               i, 
               (*self)->weightDstCellIds[i], (*self)->weightDstFaceEdgeIds[i], 
               (*self)->weightSrcCellIds[i], (*self)->weightSrcFaceEdgeIds[i],
               (*self)->weights[i]);
    }
    return 0;
}

