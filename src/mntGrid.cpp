
#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>
#include <sstream>
#include <limits>

#include "mntLogger.h"
#include <mntGrid.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGridWriter.h>
#include <fstream>
#include <netcdf.h>
#include <string>
#include <cstring>
#include <algorithm>
#include <array>
#include "mntUgrid2D.h"
#include "mntFileMeshNameExtractor.h"


/**
 * Get the paralleleliped 2D area
 * @param p0 base point
 * @param p1 point 1
 * @param p2 point 2
 * @return area
 */
inline double getArea2D(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
    Vec3 d10 = p1 - p0;
    Vec3 d20 = p2 - p0;
    return d10[0]*d20[1] - d10[1]*d20[0];
}
inline double getArea2D(const double* p0, const double* p1, const double* p2) {
    return (p1[0]-p0[0])*(p2[1] - p0[1]) - (p1[1] - p0[1])*(p2[0] - p0[0]);
}

/**
 * Fix the longitude by adding/subtracting a period to reduce the edge lengths
 * @param periodX periodicity length in x
 * @param lonBase base/reference longitude
 * @param lon longitude
 * @return corrected longitude
 */
inline
double fixLongitude(double periodX, double lonBase, double lon) {

    double diffLon = lon - lonBase;
    const double eps = 100 * std::numeric_limits<double>::epsilon();

    // favour leaving lon as is if lon == 0
    std::vector<double> diffLonMinusZeroPlus{std::abs(diffLon - periodX - eps),
                                             std::abs(diffLon),
                                             std::abs(diffLon + periodX + eps)};

    std::vector<double>::iterator it = std::min_element(diffLonMinusZeroPlus.begin(), diffLonMinusZeroPlus.end());
    int indexMin = (int) std::distance(diffLonMinusZeroPlus.begin(), it);

    // fix the longitude
    return lon + (indexMin - 1)*periodX;
}


/**
 *  Get the cell points, without applying any regularization
 * @param cellId cell ID
 * @param xyz coordinates points
 * @param face2nodes face (cell) to nodes (points) connectivity
 * @param nodes array of coordinates (output)
 */
void getCellPoints(std::size_t cellId,
                   const double xyz[], 
                   const std::size_t face2nodes[], 
                   std::vector<Vec3>& nodes) {

    const std::size_t* ptIds = &face2nodes[cellId*MNT_NUM_VERTS_PER_QUAD];

    for (std::size_t i = 0; i < MNT_NUM_VERTS_PER_QUAD; ++i) {

        std::size_t ptId = ptIds[i];
        for (auto j = 0; j < 3; ++j) {
            nodes[i][j] = xyz[ptId*3 + j];
        }
    }
}


/**
 *  Compute the face (cell) to edges connectivity
 * @param xyz coordinates points
 * @param faceNodeConnectivity face (cell) to nodes (points) connectivity
 * @param edgeNodeConnectivity edge to nodes (points) connectivity
 * @param faveEdgeConnectivity (output)
 */
void computeFaceEdgeConnectivity(const std::vector<std::size_t>& faceNodeConnectivity, 
                                 const std::vector<std::size_t>& edgeNodeConnectivity,
                                 std::vector<std::size_t>& faceEdgeConnectivity) {

    std::size_t ncells = faceNodeConnectivity.size() / MNT_NUM_VERTS_PER_QUAD;
    std::size_t nedges = edgeNodeConnectivity.size() / MNT_NUM_VERTS_PER_EDGE;

    std::map< std::array<std::size_t, 2>, std::size_t > node2Edge;
    for (std::size_t iedge = 0; iedge < nedges; ++iedge) {
        // start node
        std::size_t n0 = edgeNodeConnectivity[iedge*2 + 0];
        // end node
        std::size_t n1 = edgeNodeConnectivity[iedge*2 + 1];
        // create two entries n0 -> n1 and n1 -> n0
        std::pair< std::array<std::size_t, 2>, std::size_t > ne1({n0, n1}, iedge);
        std::pair< std::array<std::size_t, 2>, std::size_t > ne2({n1, n0}, iedge);
        node2Edge.insert(ne1);
        node2Edge.insert(ne2);
    }

    faceEdgeConnectivity.resize(ncells * MNT_NUM_EDGES_PER_QUAD);
    for (std::size_t icell = 0; icell < ncells; ++icell) {
        for (std::size_t i0 = 0; i0 < MNT_NUM_VERTS_PER_QUAD; ++i0) {
            std::size_t i1 = (i0 + 1) % MNT_NUM_VERTS_PER_QUAD;
            // start and end node indices
            std::size_t n0 = faceNodeConnectivity[icell*MNT_NUM_VERTS_PER_QUAD + i0];
            std::size_t n1 = faceNodeConnectivity[icell*MNT_NUM_VERTS_PER_QUAD + i1];
            std::size_t edgeId = node2Edge[std::array<std::size_t, 2>{n0, n1}];
            // set the edge Id for these two nodes
            faceEdgeConnectivity[icell*MNT_NUM_EDGES_PER_QUAD + i0] = edgeId;
        }
    }
}


LIBRARY_API
int mnt_grid_new(Grid_t** self) {

    *self = new Grid_t();
    (*self)->pointData = NULL;
    (*self)->points = NULL;
    (*self)->grid = NULL;
    (*self)->reader = NULL;
    (*self)->doubleArrays.resize(0);
    (*self)->numEdges = 0;

    (*self)->fixLonAcrossDateline = true;
    (*self)->averageLonAtPole = true;
    (*self)->periodX = 360.0; // if in radians, only used if the above switches are set
    (*self)->verts = NULL;
    (*self)->ownsVerts = false;
    (*self)->degrees = false;

    return 0;
}

LIBRARY_API
int mnt_grid_del(Grid_t** self) {

    for (std::size_t i = 0; i < (*self)->doubleArrays.size(); ++i) {
        (*self)->doubleArrays[i]->Delete();
    }
    if ((*self)->reader) {
        // the grid was created by reding from a file
        (*self)->reader->Delete();
    }
    else {
        if ((*self)->pointData) (*self)->pointData->Delete();
        if ((*self)->points) (*self)->points->Delete();
        if ((*self)->grid) (*self)->grid->Delete();
    }

    if ((*self)->ownsVerts) {
        delete[] (*self)->verts;
    }

    delete *self;

    return 0;
}

LIBRARY_API
int mnt_grid_setFlags(Grid_t** self, int fixLonAcrossDateline, int averageLonAtPole, int degrees) {

    (*self)->fixLonAcrossDateline = true;
    if (fixLonAcrossDateline == 0) {
        (*self)->fixLonAcrossDateline = false;
    }

    (*self)->averageLonAtPole = true;
    if (averageLonAtPole == 0) {
        (*self)->averageLonAtPole = false;
    }

    (*self)->periodX = 360;
    if (degrees == 0) {
      // lon-lat are in radians
      (*self)->periodX = 2 * M_PI;
    }

    (*self)->degrees = false;
    if (degrees != 0) {
        (*self)->degrees = true;
    }

    return 0;
}

LIBRARY_API
int mnt_grid_setPointsPtr(Grid_t** self, double points[]) {

    (*self)->verts = points;
    (*self)->ownsVerts = false;
    return 0;
}

LIBRARY_API
int mnt_grid_getPointsPtr(Grid_t** self, double** pointsPtr) {

    if ((*self)->verts) {
        *pointsPtr = (*self)->verts;
    }
    else {
        return 1; // the points have not been set
    }
    return 0;
}


LIBRARY_API
int mnt_grid_build(Grid_t** self, int nVertsPerCell, vtkIdType ncells) {

    std::string msg;

    if ((*self)->grid && (*self)->points && (*self)->pointData) {
        msg = "VTK grid has already been built, bailing out";
        mntlog::info(__FILE__, __func__, __LINE__, msg);
        return 0;
    }

    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();

    if (!(*self)->verts) {
        msg = "must call setPointsPtr before invoking build";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 2;
    }

    int save = 1;
    int npoints = nVertsPerCell * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetNumberOfComponents(3);
    // this must be called after setPointsPtr
    (*self)->pointData->SetVoidArray((*self)->verts, npoints*3, save);

    (*self)->points->SetData((*self)->pointData);

    (*self)->grid->Initialize();
    (*self)->grid->AllocateExact(ncells, MNT_NUM_VERTS_PER_QUAD);

    int cellType = -1;
    if (nVertsPerCell == MNT_NUM_VERTS_PER_QUAD) {
        cellType = VTK_QUAD;
    }
    else if (nVertsPerCell == 8) {
        cellType = VTK_HEXAHEDRON;
    }
    else {
        // error
        std::string msg = "unsupported cell type";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
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

    (*self)->grid->Squeeze();

    // clean
    ptIds->Delete();

    return 0;
}

LIBRARY_API
int mnt_grid_attach(Grid_t** self, const char* varname, int nDataPerCell, const double data[]) {

    if (!(*self)->grid) {
        std::string msg = "must build cell by cell grid before attaching a field";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
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

LIBRARY_API
int mnt_grid_computeEdgeArcLengths(Grid_t** self) {

    if (!(*self)->grid) {
        std::string msg = "must build cell by cell grid before computing arc lengths";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    std::size_t numCells;
    mnt_grid_getNumberOfCells(self, &numCells);

    if ((*self)->edgeArcLengths.size() == numCells * MNT_NUM_EDGES_PER_QUAD) {
        // already done
        return 0;
    }

    Vec3 p0, p1;

    (*self)->edgeArcLengths.resize(numCells * MNT_NUM_EDGES_PER_QUAD);
    for (std::size_t cellId = 0; cellId < numCells; ++cellId) {
        for (int edgeIndex = 0; edgeIndex < MNT_NUM_EDGES_PER_QUAD; ++edgeIndex) {

            vtkIdType ptId0 = MNT_NUM_VERTS_PER_QUAD*cellId + edgeIndex;
            vtkIdType ptId1 = MNT_NUM_VERTS_PER_QUAD*cellId + (edgeIndex + 1) % MNT_NUM_VERTS_PER_QUAD;

            (*self)->points->GetPoint(ptId0, &p0[0]);
            (*self)->points->GetPoint(ptId1, &p1[0]);

            // assumes points are in degrees
            double lam0 = p0[LON_INDEX] * M_PI/180.;
            double the0 = p0[LAT_INDEX] * M_PI/180.;
            double lam1 = p1[LON_INDEX] * M_PI/180.;
            double the1 = p1[LAT_INDEX] * M_PI/180.;

            double cos_lam0 = cos(lam0);
            double sin_lam0 = sin(lam0);
            double cos_the0 = cos(the0);
            double sin_the0 = sin(the0);

            double cos_lam1 = cos(lam1);
            double sin_lam1 = sin(lam1);
            double cos_the1 = cos(the1);
            double sin_the1 = sin(the1);

            // edge length is angle between the two points. Assume radius = 1. Angle is
            // acos of dot product in Cartesian space.
            double r0DotR1 = cos_the0*cos_lam0*cos_the1*cos_lam1 +
                             cos_the0*sin_lam0*cos_the1*sin_lam1 +
                             sin_the0*sin_the1;

            std::size_t k = MNT_NUM_EDGES_PER_QUAD*cellId + edgeIndex;
            (*self)->edgeArcLengths[k] = std::abs( acos(r0DotR1) );
        }
    }

    // add the field
    int ier = mnt_grid_attach(self, (*self)->EDGE_LENGTH_NAME.c_str(), 
                              MNT_NUM_EDGES_PER_QUAD, &(*self)->edgeArcLengths[0]);

    return ier;
}

LIBRARY_API
int mnt_grid_getEdgeArcLength(Grid_t** self, vtkIdType cellId, int edgeIndex, double* res) {

    if ((*self)->edgeArcLengths.size() == 0) {
        std::string msg = 
            "need to call mnt_grid_computeEdgeArcLengths before invoking getEdgeArcLength";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }
    std::size_t k = cellId*MNT_NUM_EDGES_PER_QUAD + edgeIndex;
    *res = (*self)->edgeArcLengths[k];
    return 0;
}

LIBRARY_API
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr) {

    if (!(*self)->grid) {
        std::string msg = "the cell by cell grid was not built";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    *grid_ptr = (*self)->grid;
    return 0;
}

LIBRARY_API
int mnt_grid_loadFromUgrid2DData(Grid_t** self, std::size_t ncells, std::size_t nedges, std::size_t npoints, 
                                 const double xyz[], const std::size_t face2nodes[], const std::size_t edge2nodes[]) {

    if ((*self)->numEdges > 0) {
        // data have already been loaded
        std::string msg = "data have already been loaded from ugrid data and the grid has been built";
        mntlog::info(__FILE__, __func__, __LINE__, msg);
        return 0;
    }

    (*self)->numEdges = nedges;

    // copy the connectivity arrays and coordinates
    (*self)->faceNodeConnectivity.resize(ncells*MNT_NUM_VERTS_PER_QUAD);
    for (std::size_t i = 0; i < ncells*MNT_NUM_VERTS_PER_QUAD; ++i) {
        (*self)->faceNodeConnectivity[i] = face2nodes[i];
    }

    (*self)->edgeNodeConnectivity.resize(nedges*MNT_NUM_VERTS_PER_EDGE);
    for (std::size_t  i = 0; i < nedges*MNT_NUM_VERTS_PER_EDGE; ++i) {
        (*self)->edgeNodeConnectivity[i] = edge2nodes[i];
    }

    // compute the face to edge connectivity from the edge-node and face-node connectivity
    computeFaceEdgeConnectivity((*self)->faceNodeConnectivity, 
                                (*self)->edgeNodeConnectivity,
                                (*self)->faceEdgeConnectivity);

    // repackage the cell vertices as a flat array
    if (npoints > 0 && (*self)->faceNodeConnectivity.size() > 0) {

        std::vector<Vec3> nodes(MNT_NUM_VERTS_PER_QUAD);

        // allocate the vertices and set the values
        (*self)->verts = new double[ncells * MNT_NUM_VERTS_PER_QUAD * 3];
        (*self)->ownsVerts = true;
        double lonBase;

        for (std::size_t icell = 0; icell < ncells; ++icell) {

            getCellPoints(icell, xyz, face2nodes, nodes);

            // fix longitude when crossing the dateline
            // use the first longitude as the base
            lonBase = nodes[0][LON_INDEX];
            // lonBase = 0;

            double avgLon = 0;
            long long poleNodeIdx = -1;
            int count = 0;
            for (auto nodeIdx = 0; nodeIdx < MNT_NUM_VERTS_PER_QUAD; ++nodeIdx) {

                double lon = nodes[nodeIdx][LON_INDEX];
                double lat = nodes[nodeIdx][LAT_INDEX];

                // lonBase *= double(nodeIdx);
                // lonBase += lon;
                // lonBase /= double(nodeIdx + 1);

                if ((*self)->fixLonAcrossDateline) {
                    lon = fixLongitude((*self)->periodX, lonBase, lon);
                }

                const double eps = 100 * std::numeric_limits<double>::epsilon();
                if (std::fabs(std::abs(lat) - 0.25*(*self)->periodX) < eps) {
                    // at the pole
                    poleNodeIdx  = nodeIdx;
                }
                else {
                    avgLon += lon;
                    count++;
                }

                // even in 2d we have three components
                (*self)->verts[LON_INDEX + nodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] = lon;
                (*self)->verts[LAT_INDEX + nodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] = lat;
                (*self)->verts[ELV_INDEX + nodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] = 0.0;
            }
            avgLon /= count;

            // check if one of the nodes is at the north/south pole. In
            // this case the longitude is ill-defined. Set the longitude there to be the
            // average of the 3 other longitudes.

            if ((*self)->averageLonAtPole && poleNodeIdx >= 0) {
                std::stringstream msg;
                msg << "setting longitude = " << (*self)->verts[LON_INDEX + poleNodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] 
                    << " in cell " << icell << " to " << avgLon << ", other lons = ";
                    for (auto j = 0; j < MNT_NUM_VERTS_PER_QUAD; ++j) {
                        msg << (*self)->verts[LON_INDEX + j*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] << ", ";
                    }
                    msg << '\n';
                mntlog::info(__FILE__, __func__, __LINE__, msg.str());
                (*self)->verts[LON_INDEX + poleNodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] = avgLon;
            }

            // run through another date line correction
            if ((*self)->fixLonAcrossDateline) {
                lonBase = (*self)->verts[LON_INDEX + 0*3 + icell*MNT_NUM_VERTS_PER_QUAD*3];
                for (auto nodeIdx = 1; nodeIdx < MNT_NUM_VERTS_PER_QUAD; ++nodeIdx) {
                    double lon = (*self)->verts[LON_INDEX + nodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3];
                    lon = fixLongitude((*self)->periodX, lonBase, lon);
                    (*self)->verts[LON_INDEX + nodeIdx*3 + icell*MNT_NUM_VERTS_PER_QUAD*3] = lon;
                }
            }
        }
    }

    // build the connectivity
    int ier = mnt_grid_build(self, MNT_NUM_VERTS_PER_QUAD, ncells);

    return ier;
}


LIBRARY_API
int mnt_grid_loadFromUgrid2DFile(Grid_t** self, const char* fileAndMeshName) {

    // extract the filename and the mesh name from "filename$meshname"
    auto fm = fileMeshNameExtractor(fileAndMeshName);

    std::string filename = fm.first;
    std::string meshname = fm.second;

    Ugrid2D ugrid;
    int ier = ugrid.load(filename, meshname);
    if (ier != 0) {
        mntlog::error(__FILE__, __func__, __LINE__, 
            "could not read mesh \"" + meshname + "\" in UGRID file \"" + filename + "\"");
        return 1;
    }

    std::size_t ncells = ugrid.getNumberOfFaces();
    std::size_t nedges = ugrid.getNumberOfEdges();
    std::size_t npoints = ugrid.getNumberOfPoints();

    const std::vector<double>& xyz = ugrid.getPoints();
    const std::vector<std::size_t>& face2nodes = ugrid.getFacePointIds();
    const std::vector<std::size_t>& edge2nodes = ugrid.getEdgePointIds();

    ier = mnt_grid_loadFromUgrid2DData(self, ncells, nedges, npoints, 
                                       &xyz[0], &face2nodes[0], &edge2nodes[0]);

    return ier;
}

LIBRARY_API
int mnt_grid_load(Grid_t** self, const char* filename) {

    // check if the file exists
    if (!fstream(filename).good()) {
        mntlog::error(__FILE__, __func__, __LINE__, 
                     "file " + std::string(filename) + "does not exist");
        return 1;
    }

    if ((*self)->grid) {
        (*self)->grid->Delete();
    }

    (*self)->reader = vtkUnstructuredGridReader::New();
    (*self)->reader->SetFileName(filename);
    (*self)->reader->Update();
    (*self)->grid = (*self)->reader->GetOutput();
    (*self)->points = (*self)->grid->GetPoints();
    (*self)->verts = (double *) (*self)->points->GetData()->GetVoidPointer(0);
    return 0;
}

LIBRARY_API
int mnt_grid_dump(Grid_t** self, const char* filename) {

    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
#if (VTK_MAJOR_VERSION >= 9)
    writer->SetFileVersion(42); // want to write old legacy format
#endif
    writer->SetFileName(filename);
    writer->SetInputData((*self)->grid);
    writer->Update();
    writer->Delete();
    return 0;
}

LIBRARY_API
int mnt_grid_print(Grid_t** self) {

    vtkPoints* points = (*self)->grid->GetPoints();
    vtkIdType npoints = points->GetNumberOfPoints();
    std::cout << "Number of points: " << npoints << '\n';

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();
    std::cout << "Number of cells: " << ncells << '\n';

    Vec3 pt;

    for (vtkIdType i = 0; i < ncells; ++i) {

        vtkCell* cell = (*self)->grid->GetCell(i);

        for (int j = 0; j < cell->GetNumberOfPoints(); ++j) {
            vtkIdType k = cell->GetPointId(j);
            (*self)->points->GetPoint(k, &pt[0]);
            std::cout << "\t point " << pt[0] << ',' << pt[1] << ',' << pt[2] << std::endl;
        }
    }

    for (std::size_t i = 0; i < (*self)->faceNodeConnectivity.size(); 
            i += MNT_NUM_VERTS_PER_QUAD) {
        auto icell = i / MNT_NUM_VERTS_PER_QUAD;
        std::cout << "\t face " << icell << ": nodes";
        for (auto k = 0; k < MNT_NUM_VERTS_PER_QUAD; ++k) {
            std::cout << " " << (*self)->faceNodeConnectivity[icell*MNT_NUM_VERTS_PER_QUAD + k];
        }
        std::cout << std::endl;
    }
    for (std::size_t i = 0; i < (*self)->edgeNodeConnectivity.size();
             i += MNT_NUM_VERTS_PER_EDGE) {
        auto iedge = i / MNT_NUM_VERTS_PER_EDGE;
        std::cout << "\t edge " << iedge << ": nodes";
        for (auto k = 0; k < MNT_NUM_VERTS_PER_EDGE; ++k) {
            std::cout << " " << (*self)->edgeNodeConnectivity[iedge*MNT_NUM_VERTS_PER_EDGE + k];
        }
        std::cout << std::endl;
    }

    return 0;
}

LIBRARY_API
int mnt_grid_getPoints(Grid_t** self, vtkIdType cellId, int edgeIndex,
                       double point0[], double point1[]) {

    // flat index for the start point, 3d coordinates
    std::size_t k0 = MNT_NUM_VERTS_PER_QUAD*3*cellId + 3*((edgeIndex + 0) % MNT_NUM_VERTS_PER_QUAD);

    // flat index for the end point, 3d coordinates
    std::size_t k1 = MNT_NUM_VERTS_PER_QUAD*3*cellId + 3*((edgeIndex + 1) % MNT_NUM_VERTS_PER_QUAD);

    if (edgeIndex < 2) {
        // edge's direction is counterclockwise
        for (std::size_t i = 0; i < 3; ++i) {
            point0[i] = (*self)->verts[i + k0];
            point1[i] = (*self)->verts[i + k1];
        }
    }
    else {
        // edge's direction is clockwise - reverse order of point0 and point1
        for (std::size_t i = 0; i < 3; ++i) {
            point1[i] = (*self)->verts[i + k0];
            point0[i] = (*self)->verts[i + k1];
        }
    }

    return 0;
}

LIBRARY_API
int mnt_grid_getNodeIds(Grid_t** self, vtkIdType cellId, int edgeIndex, vtkIdType nodeIds[]) {

    if ((*self)->faceNodeConnectivity.size() == 0) {
        std::string msg = "no face-node connectivity, grid is empty or was not built from Ugrid";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    // nodeIndex0,1 are the local cell indices of the vertices in the range 0-3
    int nodeIndex0 = edgeIndex;
    int nodeIndex1 = (edgeIndex + 1) % MNT_NUM_VERTS_PER_QUAD;

    // edges 2-3 go clockwise
    // edges 0-1 go anticlockwise
    if (edgeIndex >= 2) {
        // swap order
        int tmp = nodeIndex0;
        nodeIndex0 = nodeIndex1;
        nodeIndex1 = tmp;
    }

    nodeIds[0] = (*self)->faceNodeConnectivity[MNT_NUM_VERTS_PER_QUAD*cellId + nodeIndex0];
    nodeIds[1] = (*self)->faceNodeConnectivity[MNT_NUM_VERTS_PER_QUAD*cellId + nodeIndex1];

    return 0;
}

LIBRARY_API
int mnt_grid_getEdgeId(Grid_t** self, vtkIdType cellId, int edgeIndex, 
                       std::size_t* edgeId, int* signEdge) {

    if ((*self)->faceNodeConnectivity.size() == 0) {
        std::string msg = "no face-node connectivity, grid is empty or was not built from Ugrid";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    // initialize
    *signEdge = 0;

    // fetch the node Ids of this edge
    vtkIdType nodeIds[2];
    int ier = mnt_grid_getNodeIds(self, cellId, edgeIndex, nodeIds);
    if (ier != 0) {
        return 1;
    }

    // iterate over the edges of this face until we find the edge
    // that has vertices nodeIds (but not necessarily in the same
    // order)
    for (int ie = 0; ie < MNT_NUM_EDGES_PER_QUAD; ++ie) {

        // edgeId under consideration
        vtkIdType eId = (*self)->faceEdgeConnectivity[MNT_NUM_EDGES_PER_QUAD*cellId + ie];

        // vertex Ids of the edge
        vtkIdType nId0 = (*self)->edgeNodeConnectivity[eId*2 + 0];
        vtkIdType nId1 = (*self)->edgeNodeConnectivity[eId*2 + 1];

        if (nId0 == nodeIds[0] && nId1 == nodeIds[1]) {
            // found edge and the direction is left->right, bottom -> up
            *signEdge = 1;
            *edgeId = eId;
            break;
        }
        else if (nId0 == nodeIds[1] && nId1 == nodeIds[0]) {
            // found edge and the direction is opposite
            *signEdge = -1;
            *edgeId = eId;
            break;
        }
    }

    return 0;
}


LIBRARY_API
int mnt_grid_getNumberOfCells(Grid_t** self, std::size_t* numCells) {

    *numCells = (*self)->grid->GetNumberOfCells();
    return 0;
}

LIBRARY_API
int mnt_grid_getNumberOfEdges(Grid_t** self, std::size_t* numEdges) {

    if ((*self)->faceNodeConnectivity.size() == 0) {
        std::string msg = "no face-node connectivity, grid is empty or was not built from Ugrid";
        mntlog::warn(__FILE__, __func__, __LINE__, msg);
    }

    *numEdges = (*self)->numEdges;
    return 0;
}

LIBRARY_API
int mnt_grid_check(Grid_t** self, std::size_t* numBadCells) {

    *numBadCells = 0;

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();

    std::vector<Vec3> pts(MNT_NUM_VERTS_PER_QUAD); // 2D

    for (vtkIdType i = 0; i < ncells; ++i) {

        vtkCell* cell = (*self)->grid->GetCell(i);

        for (int j = 0; j < cell->GetNumberOfPoints(); ++j) {
            vtkIdType k = cell->GetPointId(j);
            (*self)->points->GetPoint(k, &pts[j][0]);
        }

        std::size_t badCell = 0;
        double area = getArea2D(pts[0], pts[1], pts[3]);
        if (area < 0.) {
            std::stringstream ss;
            ss << pts[0] << ';' << pts[1] << ';' << pts[3];
            mntlog::warn(__FILE__, __func__, __LINE__, 
            "cell " + std::to_string(i) + 
            " has negative area = " + std::to_string(area) +
            " for points 0-1-3 " + ss.str() + "\n");
            badCell = 1;            
        }

        area = getArea2D(pts[2], pts[3], pts[1]);
        if (area < 0.) {
            std::stringstream ss;
            ss << pts[2] << ';' << pts[3] << ';' << pts[1];
            mntlog::warn(__FILE__, __func__, __LINE__, 
            "cell " + std::to_string(i) + 
            " has negative area = " + std::to_string(area) +
            " for points 2-3-1 " + ss.str() + "\n");
            badCell = 1;            
        }
        (*numBadCells) += badCell;
    }

    return 0;
}
