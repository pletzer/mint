#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <string>
#include <limits>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkUnstructuredGrid.h>
#include <mntVecN.h>

#ifndef MNT_UGRID_2D
#define MNT_UGRID_2D

#define NUM_PARAM_DIMS 2
#define NUM_SPACE_DIMS 3

class Ugrid2D {

public:

/**
 * Constructor
 */
Ugrid2D() {
	this->isCartesian = false;
}

/**
 * Destructor
 */
~Ugrid2D() {
}


/**
 * Load from Ugrid file 
 * @param filename file name
 * @param meshname mesh name
 * @return error (0=OK)
 */
int load(const std::string& filename, const std::string& meshname);


/**
 * Dump the grid to a Vtk file
 * @param filename file name
 */
void dumpGridVtk(const std::string& filename);

/**
 * Get the number of faces
 * @return number
 */
size_t getNumberOfFaces() const {
    return this->numFaces;
}

/**
 * Get the number of unique edges
 * @return number
 */
size_t getNumberOfEdges() const {
    return this->numEdges;
}

/**
 * Get the number of points/vertices
 * @return number
 */
size_t getNumberOfPoints() const {
    return this->numPoints;
}

/**
 * Get pointer to the point Ids
 * @param face Id
 * @return pointer
 */
const size_t* getFacePointIds(size_t faceId) const {
    return &this->face2Points[faceId*4];
}
const std::vector<size_t>& getFacePointIds() const {
    return this->face2Points;
}

/**
 * Get pointer to the edge Ids
 * @param face Id 
 * @return pointer
 */
const size_t* getFaceEdgeIds(size_t faceId) const {
    return &this->face2Edges[faceId*4];
}
const std::vector<size_t>& getFaceEdgeIds() const {
    return this->face2Edges;
}

/**
 * Get pointer to the point Ids
 * @param edge Id 
 * @return pointer
 */
const size_t* getEdgePointIds(size_t edgeId) const {
    return &this->edge2Points[edgeId*2];
}
const std::vector<size_t>& getEdgePointIds() const {
    return this->edge2Points;
}

/**
 * Get pointer to the coordinates
 * @param pointId point id
 * @return pointer
 */
const double* getPoint(size_t pointId) const {
    return &this->points[pointId*NUM_SPACE_DIMS];
}
const std::vector<double>& getPoints() const {
    return this->points;
}

/**
 * Get the regularized points of the edge
 * @param edgeId edge Id
 * @return the two points with 360 added/substracted to minimize the edge length
 */
std::vector<Vec3> getEdgePointsRegularized(size_t edgeId) const;

/**
 * Get the regularized points of the face
 * @param faceId face Id
 * @return the four points with 360 added/substracted to make the cell have a positive area
 */
std::vector<Vec3> getFacePointsRegularized(size_t faceId) const;

/**
 * Get min/max range of the domain
 * @param xmin low point of the domain (output)
 * @param xmax high point of the domain (output)
 */
void getRange(double xmin[], double xmax[]) const;


/**
 * Get the face vertex coordinates
 * @return array of points
 */
std::vector<Vec3> getFacePoints(size_t faceId) const;

/**
 * Get the edge vertex coordinates
 * @return array of points
 */
std::vector<Vec3> getEdgePoints(size_t edgeId) const;


private:

    // vertex coordinates
    std::vector<double> points;

    // domain min/max
    Vec3 xmin;
    Vec3 xmax;

    size_t numPoints;
    size_t numEdges;
    size_t numFaces;

    // face to node connectivity
    std::vector<size_t> face2Points;

    // face to edge connectivity
    std::vector<size_t> face2Edges;

    // edge to node connectivity
    std::vector<size_t> edge2Points;

    int readConnectivityData(int ncid, int meshid, 
                             const std::string& role,
                             std::vector<size_t>& data);

    int readPoints(int ncid, int meshid);

    void fixPeriodicity();

    bool isCartesian;

};

#endif // MNT_UGRID_2D
