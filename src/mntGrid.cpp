
#define _USE_MATH_DEFINES // M_PI for Visual Studio
#include <cmath>
#include <sstream>

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
 * Fix the longitude by adding/subtracting a period to reduce the edge lengths
 * @param periodX periodicity length in x
 * @param lonBase base/reference longitude
 * @param lon longitude
 * @return corrected longitude
 */
inline
double fixLongitude(double periodX, double lonBase, double lon) {

    double diffLon = lon - lonBase;

    std::vector<double> diffLonMinusZeroPlus{std::abs(diffLon - periodX),
                                             std::abs(diffLon),
                                             std::abs(diffLon + periodX)};

    std::vector<double>::iterator it = std::min_element(diffLonMinusZeroPlus.begin(), diffLonMinusZeroPlus.end());
    int indexMin = (int) std::distance(diffLonMinusZeroPlus.begin(), it);

    // fix the longitude
    return lon + (indexMin - 1)*periodX;
}

LIBRARY_API
int mnt_grid_new(Grid_t** self) {

    *self = new Grid_t();
    (*self)->pointData = NULL;
    (*self)->points = NULL;
    (*self)->grid = NULL;
    (*self)->reader = NULL;
    (*self)->doubleArrays.resize(0);

    (*self)->fixLonAcrossDateline = true;
    (*self)->averageLonAtPole = true;
    (*self)->periodX = 360.0; // if in radians, only used if the above switches are set
    (*self)->verts = NULL;
    (*self)->ownsVerts = false;

    return 0;
}

LIBRARY_API
int mnt_grid_del(Grid_t** self) {

    for (std::size_t i = 0; i < (*self)->doubleArrays.size(); ++i) {
        (*self)->doubleArrays[i]->Delete();
    }
    if ((*self)->reader) {
        (*self)->reader->Delete();
    }
    else {
        if ((*self)->grid) (*self)->grid->Delete();
    }
    if ((*self)->points) (*self)->points->Delete();
    if ((*self)->pointData) (*self)->pointData->Delete();

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

    return 0;
}

LIBRARY_API
int mnt_grid_setPointsPtr(Grid_t** self, double points[]) {

    (*self)->verts = points;
    (*self)->ownsVerts = false;
    return 0;
}

LIBRARY_API
int mnt_grid_build(Grid_t** self, int nVertsPerCell, vtkIdType ncells) {

    (*self)->pointData = vtkDoubleArray::New();
    (*self)->points = vtkPoints::New();
    (*self)->grid = vtkUnstructuredGrid::New();

    int save = 1;
    int npoints = nVertsPerCell * ncells;
    (*self)->pointData->SetNumberOfTuples(npoints);
    (*self)->pointData->SetNumberOfComponents(3);
    // this must be called after setPointsPtr
    (*self)->pointData->SetVoidArray((*self)->verts, npoints*3, save);

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
    // (*self)->grid->BuildLinks(); // DO WE NEED THIS?

    // clean
    ptIds->Delete();

    return 0;
}

LIBRARY_API
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

LIBRARY_API
int mnt_grid_computeEdgeArcLengths(Grid_t** self) {

    if (!(*self)->grid) {
        return 1;
    }

    std::size_t numCells;
    mnt_grid_getNumberOfCells(self, &numCells);

    if ((*self)->edgeArcLengths.size() == numCells * 4) {
        // already done
        return 0;
    }

    Vec3 p0, p1;

    (*self)->edgeArcLengths.resize(numCells * 4);
    for (std::size_t cellId = 0; cellId < numCells; ++cellId) {
        for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex) {

            vtkIdType ptId0 = 4*cellId + edgeIndex;
            vtkIdType ptId1 = 4*cellId + (edgeIndex + 1) % 4;

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
            double r0DotR1 = cos_the0*cos_lam0*cos_the1*cos_lam1 + cos_the0*sin_lam0*cos_the1*sin_lam1 + sin_the0*sin_the1;
            std::size_t k = 4*cellId + edgeIndex;
            (*self)->edgeArcLengths[k] = std::abs( acos(r0DotR1) );
        }
    }

    // add the field
    int ier = mnt_grid_attach(self, (*self)->EDGE_LENGTH_NAME.c_str(), 4, &(*self)->edgeArcLengths[0]);

    return ier;
}

LIBRARY_API
int mnt_grid_getEdgeArcLength(Grid_t** self, vtkIdType cellId, int edgeIndex, double* res) {

    if ((*self)->edgeArcLengths.size() == 0) {
        mntlog::error(__FILE__, __func__, __LINE__, 
            "you need to call mnt_grid_computeEdgeArcLengths before invoking mnt_grid_getEdgeArcLength");
        return 1;
    }
    std::size_t k = cellId*4 + edgeIndex;
    *res = (*self)->edgeArcLengths[k];
    return 0;
}

LIBRARY_API
int mnt_grid_get(Grid_t** self, vtkUnstructuredGrid** grid_ptr) {
    *grid_ptr = (*self)->grid;
    return 0;
}

LIBRARY_API
int mnt_grid_loadFrom2DUgrid(Grid_t** self, const char* fileAndMeshName) {

    // extract the filename and the mesh name from "filename:meshname"
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

    double xmin[3], xmax[3];
    ugrid.getRange(xmin, xmax);
    double lonMin = xmin[0];

    // copy
    (*self)->faceNodeConnectivity = ugrid.getFacePointIds();
    (*self)->edgeNodeConnectivity = ugrid.getEdgePointIds();

    std::size_t ncells = ugrid.getNumberOfFaces();
    std::size_t nedges = ugrid.getNumberOfEdges();
    std::size_t npoints = ugrid.getNumberOfPoints();
    int numVertsPerCell = 4;

    // get the face to edge connectivity from the file
    (*self)->faceEdgeConnectivity = ugrid.getFaceEdgeIds();

    if ((*self)->faceEdgeConnectivity.size() == 0) {

        // compute the face to edge connectivity from the edge-node and face-node connectivity
        std::map< std::array<std::size_t, 2>, std::size_t > node2Edge;
        for (std::size_t iedge = 0; iedge < nedges; ++iedge) {
            // start node
            std::size_t n0 = (*self)->edgeNodeConnectivity[iedge*2 + 0];
            // end node
            std::size_t n1 = (*self)->edgeNodeConnectivity[iedge*2 + 1];
            // create two entries n0 -> n1 and n1 -> n0
            std::pair< std::array<std::size_t, 2>, std::size_t > ne1({n0, n1}, iedge);
            std::pair< std::array<std::size_t, 2>, std::size_t > ne2({n1, n0}, iedge);
            node2Edge.insert(ne1);
            node2Edge.insert(ne2);
        }
        (*self)->faceEdgeConnectivity.resize(ncells * 4);
        for (std::size_t icell = 0; icell < ncells; ++icell) {
            for (std::size_t i0 = 0; i0 < 4; ++i0) {
                std::size_t i1 = (i0 + 1) % 4;
                // start and end node indices
                std::size_t n0 = (*self)->faceNodeConnectivity[icell*4 + i0];
                std::size_t n1 = (*self)->faceNodeConnectivity[icell*4 + i1];
                std::size_t edgeId = node2Edge[std::array<std::size_t, 2>{n0, n1}];
                // set the edge Id for these two nodes
                (*self)->faceEdgeConnectivity[icell*4 + i0] = edgeId;
            }
        }
    }

    // repackage the cell vertices as a flat array

    if (npoints > 0 && (*self)->faceNodeConnectivity.size() > 0) {

        // allocate the vertices and set the values
        (*self)->verts = new double[ncells * numVertsPerCell * 3];
        (*self)->ownsVerts = true;

        for (std::size_t icell = 0; icell < ncells; ++icell) {

            // fix longitude when crossing the dateline
            // use the first longitude as the base
            std::size_t kBase = (*self)->faceNodeConnectivity[icell*numVertsPerCell];
            double lonBase = ugrid.getPoint(kBase)[LON_INDEX];

            double avgLon = 0;
            long long poleNodeIdx = -1;
            int count = 0;
            for (auto nodeIdx = 0; nodeIdx < numVertsPerCell; ++nodeIdx) {

                std::size_t k = (*self)->faceNodeConnectivity[icell*numVertsPerCell + nodeIdx];
                double lon = ugrid.getPoint(k)[LON_INDEX]; //lons[k];
                double lat = ugrid.getPoint(k)[LAT_INDEX];

                if ((*self)->fixLonAcrossDateline) {
                    lon = fixLongitude((*self)->periodX, lonBase, lon);
                }

                if (std::abs(lat) == 0.25*(*self)->periodX) {
                    // at the pole
                    poleNodeIdx  = nodeIdx;
                }
                else {
                    avgLon += lon;
                    count++;
                }

                // even in 2d we have three components
                (*self)->verts[LON_INDEX + nodeIdx*3 + icell*numVertsPerCell*3] = lon;
                (*self)->verts[LAT_INDEX + nodeIdx*3 + icell*numVertsPerCell*3] = lat;
                (*self)->verts[ELV_INDEX + nodeIdx*3 + icell*numVertsPerCell*3] = 0.0;
            }
            avgLon /= count;

            // check if there if one of the cell nodes is at the north/south pole. In
            // this case the longitude is ill-defined. Set the longitude there to the
            // average of the 3 other longitudes.

            if ((*self)->averageLonAtPole && poleNodeIdx >= 0) {
                (*self)->verts[LON_INDEX + poleNodeIdx*3 + icell*numVertsPerCell*3] = avgLon;
            }

            // make sure the cell is within the lonMin to lonMin + periodX range
            double offsetLon = 0.0;
            if ((*self)->fixLonAcrossDateline) {
                if (avgLon > lonMin + (*self)->periodX) {
                    offsetLon = -(*self)->periodX;
                }
                else if (avgLon < lonMin) {
                    offsetLon = (*self)->periodX;
                }
                for (int nodeIdx = 0; nodeIdx < numVertsPerCell; ++nodeIdx) {
                    (*self)->verts[LON_INDEX + (std::size_t) nodeIdx * 3 + icell*numVertsPerCell*3] += offsetLon;
                }
            }
        }
    }

    // build the connectivity
    ier = mnt_grid_build(self, numVertsPerCell, ncells);

    return 0;
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
    return 0;
}

LIBRARY_API
int mnt_grid_dump(Grid_t** self, const char* filename) {
    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
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
            std::cout << "\tpoint " << pt[0] << ',' << pt[1] << ',' << pt[2] << '\n';
        }
    }

    return 0;
}

LIBRARY_API
int mnt_grid_getPoints(Grid_t** self, vtkIdType cellId, int edgeIndex,
                       double point0[], double point1[]) {

    // flat index for the start point, 4 points per cell, 3d coordinates
    std::size_t k0 = 4*3*cellId + 3*((edgeIndex + 0) % 4);

    // flat index for the end point, 4 points per cell, 3d coordinates
    std::size_t k1 = 4*3*cellId + 3*((edgeIndex + 1) % 4);

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

    // nodeIndex0,1 are the local cell indices of the vertices in the range 0-3
    int nodeIndex0 = edgeIndex;
    // 4 vertices per cell
    int nodeIndex1 = (edgeIndex + 1) % 4;

    // edges 2-3 go clockwise
    // edges 0-1 go anticlockwise
    if (edgeIndex >= 2) {
        // swap order
        int tmp = nodeIndex0;
        nodeIndex0 = nodeIndex1;
        nodeIndex1 = tmp;
    }

    nodeIds[0] = (*self)->faceNodeConnectivity[4*cellId + nodeIndex0];
    nodeIds[1] = (*self)->faceNodeConnectivity[4*cellId + nodeIndex1];

    return 0;
}

LIBRARY_API
int mnt_grid_getEdgeId(Grid_t** self, vtkIdType cellId, int edgeIndex, 
                       std::size_t* edgeId, int* signEdge) {

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
    for (int ie = 0; ie < 4; ++ie) {

        // edgeId under consideration
        vtkIdType eId = (*self)->faceEdgeConnectivity[4*cellId + ie];

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

    *numEdges = (*self)->edgeNodeConnectivity.size() / 2;
    return 0;
}

double getArea2D(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
    Vec3 d10 = p1 - p0;
    Vec3 d20 = p2 - p0;
    return d10[0]*d20[1] - d10[1]*d20[0];
}

LIBRARY_API
int mnt_grid_check(Grid_t** self, std::size_t* numBadCells) {

    *numBadCells = 0;

    vtkIdType ncells = (*self)->grid->GetNumberOfCells();

    std::vector<Vec3> pts(4); // 2D

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
