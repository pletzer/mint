#include <mntUgrid2D.h>
#include <mntLineLineIntersector.h>
#include <netcdf.h>
#include <iostream>
#include <set>
#include <cmath>
#include <algorithm>
#include <limits>
#include <fstream>

#define LON_INDEX 0
#define LAT_INDEX 1
#define ELV_INDEX 2

#define X_INDEX 0
#define Y_INDEX 1
#define Z_INDEX 2


struct LambdaBegFunctor {
    // compare two elements of the array
    bool operator()(const std::pair<size_t, std::vector<double> >& x, 
                    const std::pair<size_t, std::vector<double> >& y) {
        return (x.second[0] < y.second[0]);
    }
};


std::vector<Vec3> 
Ugrid2D::getFacePoints(size_t faceId) const {

    const size_t* pointIds = this->getFacePointIds(faceId);
    std::vector<Vec3> res(4); // 2d 4 points per quad

    // iterate over the 4 points
    for (size_t i = 0; i < 4; ++i) {
        const double* p = this->getPoint(pointIds[i]);
        res[i] = Vec3(p);
    }

    return res;
}


std::vector<Vec3> 
Ugrid2D::getEdgePoints(size_t edgeId) const {

    std::vector<Vec3> res;

    // itereate over the 2 points spanning the edge
    for (size_t i = 0; i < 2; ++i) {

        // get the point id
        size_t pointId = this->edge2Points[i + edgeId*2];

        // get the coordinates of this point
        const double* p = this->getPoint(pointId);

        // add
        res.push_back(Vec3(p));
    }
    return res;
}

void
Ugrid2D::getRange(double xMin[], double xMax[]) const {
    for (size_t i = 0; i < NUM_SPACE_DIMS; ++i) {
        xMin[i] = this->xmin[i];
        xMax[i] = this->xmax[i];
    }
}

int 
Ugrid2D::load(const std::string& filename, const std::string& meshname) {

    int ier = 0;
    int ncid;

    // open the file
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot open \"" << filename << "\"\n";
        return 1;
    }

    // mesh variable
    int meshid;
    ier = nc_inq_varid(ncid, meshname.c_str(), &meshid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot find variable named \"" << meshname << "\"\n";
        nc_close(ncid);
        return 2;
    }

    ier = this->readConnectivityData(ncid, meshid, 
                "face_node_connectivity", this->face2Points);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: variable \"" << meshname
        << "\" does not have attribute \"face_node_connectivity\"\n";
        nc_close(ncid);
        return 3;
    }

    ier = this->readConnectivityData(ncid, meshid, 
                "face_edge_connectivity", this->face2Edges);
    /* no longer an error
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: variable \"" << meshname
        << "\" does not have attribute \"face_edge_connectivity\"\n";
        nc_close(ncid);
        return 4;
    }
    */

    ier = this->readConnectivityData(ncid, meshid, 
                "edge_node_connectivity", this->edge2Points);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: variable \"" << meshname
        << "\" does not have attribute \"edge_node_connectivity\"\n";
        nc_close(ncid);
        return 5;
    }

    // read the node coordinates
    ier = this->readPoints(ncid, meshid);
    if (ier != NC_NOERR) {
        std::cerr << "ERROR: cannot read node coordinates for mesh \"" 
        << meshname << "\"\n";
        nc_close(ncid);
        return 6;
    }

    size_t n = this->face2Points.size();
    this->numFaces = n / 4; // 4 points per face
    n = this->edge2Points.size();
    this->numEdges = n / 2; // 2 points per edge


    // compute min/max values after regularizing the coords across the faces
    // (ie adding/subtracting 360 deg for the longitude to make the face area positive)
    this->xmin = +std::numeric_limits<double>::max();
    this->xmax = -std::numeric_limits<double>::max();
    for (size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        std::vector<Vec3> nodes = this->getFacePointsRegularized(faceId);
        for (const Vec3& p : nodes) {
            for (size_t j = 0; j < this->xmin.size(); ++j) {
                this->xmin[j] = (p[j] < this->xmin[j]? p[j]: this->xmin[j]);
                this->xmax[j] = (p[j] > this->xmax[j]? p[j]: this->xmax[j]);
            }                     
        }
    }

    return 0;
}

int
Ugrid2D::readConnectivityData(int ncid, int meshid, 
                              const std::string& role,
                              std::vector<size_t>& data) {

    int ier;

    // get the lengths of the attribute string
    size_t len;
    ier = nc_inq_attlen(ncid, meshid, role.c_str(), &len);
    if (ier != NC_NOERR) return 1;

    // read the attribute value, ie the name of the variables we will need to read
    std::string varname(len, ' ');
    ier = nc_get_att_text(ncid, meshid, role.c_str(), &varname[0]);
    if (ier != NC_NOERR) return 2;

    // fetch the variable Id for this variable name
    int varid;
    ier = nc_inq_varid(ncid, varname.c_str(), &varid);
    if (ier != NC_NOERR) return 3;

    // dimensions of the variable to read
    int dimids[2];
    size_t n0, n1;
    // either 0 or 1
    int startIndex;

    // read the data
    ier = nc_inq_vardimid(ncid, varid, dimids);
    if (ier != NC_NOERR) return 4;

    ier = nc_inq_dimlen(ncid, dimids[0], &n0);
    if (ier != NC_NOERR) return 5;

    ier = nc_inq_dimlen(ncid, dimids[1], &n1);
    if (ier != NC_NOERR) return 6;

    ier = nc_get_att_int(ncid, varid, "start_index", &startIndex);
    if (ier != NC_NOERR) return 8;

    std::vector<unsigned long long> buffer(n0 * n1);
    ier = nc_get_var_ulonglong(ncid, varid, &buffer[0]);
    if (ier != NC_NOERR) return 7;

    // subtract start_index
    data.resize(n0 * n1);
    for (size_t i = 0; i < n0 * n1; ++i) {
        data[i] = buffer[i] - startIndex;
    }

    return 0;
}

int
Ugrid2D::readPoints(int ncid, int meshid) {

    int ier;

    // get the lengths of the attribute string
    size_t len;
    ier = nc_inq_attlen(ncid, meshid, "node_coordinates", &len);
    if (ier != NC_NOERR) return 10;

    // read the attribute, lists the name of the lon and lat coordinates
    std::string val(len, ' ');
    ier = nc_get_att_text(ncid, meshid, "node_coordinates", &val[0]);
    if (ier != NC_NOERR) return 11;

    // val is "varx vary" where var{x,y} are the variable names
    size_t n = val.size();
    size_t spaceL = val.find(' ');
    size_t spaceR = val.rfind(' ');
    if (spaceL >= n) {
        // could not find space
        std::cerr << "ERROR: node_coordinates attribute \""
        << val << "\" should contain space separated list of coordinate names\n";
        return 19;
    }
    std::string varx = val.substr(0, spaceL);
    std::string vary = val.substr(spaceR + 1, n - spaceR - 1);

    int varxid, varyid;
    ier = nc_inq_varid(ncid, varx.c_str(), &varxid);
    if (ier != NC_NOERR) return 12;
    ier = nc_inq_varid(ncid, vary.c_str(), &varyid);
    if (ier != NC_NOERR) return 13;

    int dimids[1];

    this->numPoints = 0;
    this->points.resize(0);
    int varids[] = {varxid, varyid};
    for (int ivar = 0; ivar < 2; ++ivar) {

        int varid = varids[ivar];

        // get the attribute length
        ier = nc_inq_attlen(ncid, varid, "standard_name", &len);
        if (ier != NC_NOERR) {
            std::cerr << "ERROR: variable with varid = " << varid
            << " has no attribute \"standard_name\"\n";
            return 14;
        }
        std::string var_stdn(len, ' ');

        // read the attribute
        ier = nc_get_att_text(ncid, varid, "standard_name", &var_stdn[0]);
        if (ier != NC_NOERR) return 15;

        // get the dimension
        ier = nc_inq_vardimid(ncid, varid, dimids);
        if (ier != NC_NOERR) return 16;

        ier = nc_inq_dimlen(ncid, dimids[0], &this->numPoints);
        if (ier != NC_NOERR) return 17;

        // allocate/resize
        std::vector<double> data(this->numPoints);
        if (this->points.size() == 0) {
            this->points.resize(this->numPoints * NUM_SPACE_DIMS, 0.0);
        }

        // read the data
        ier = nc_get_var_double(ncid, varid, &data[0]);
        if (ier != NC_NOERR) {
            std::cerr << "ERROR: could not read \""
            << varx << "\"\n";
            return 18;
        }
        
        // associate data  our coordinate variable
        size_t j;
        if (var_stdn == "longitude") {
            j = LON_INDEX;
        }
        else if (var_stdn == "latitude") {
            j = LAT_INDEX;
        }
        else if (var_stdn == "x-coordinate in Cartesian system") {
            j = X_INDEX;
            this->isCartesian = true;
        }
        else if (var_stdn == "y-coordinate in Cartesian system") {
            j = Y_INDEX;
        }
        else if (var_stdn == "z-coordinate in Cartesian system") {
            j = Z_INDEX;
        }
        else {
            std::cerr << "ERROR: unknown coordinate with standard_name \""
            << var_stdn << "\"\n";
            return 19; 
        }
        for (size_t i = 0; i < this->numPoints; ++i) {
            this->points[j + NUM_SPACE_DIMS*i] = data[i];
        }
    }

    return 0;
}

std::vector<Vec3> 
Ugrid2D::getFacePointsRegularized(size_t faceId) const {

    std::vector<Vec3> res = this->getFacePoints(faceId);

    if (this->isCartesian) {
        return res;
    }

    // regularize
    for (size_t i = 1; i < res.size(); ++i) {

        // add/subtract 360 
        double dLon = res[i][LON_INDEX] - res[0][LON_INDEX];
        double dLonsPM360[] = {std::abs(dLon - 360.), 
                               std::abs(dLon       ), 
                               std::abs(dLon + 360.)};
        double* minDLon = std::min_element(&dLonsPM360[0], &dLonsPM360[3]);
        int indexMin = (int) std::distance(dLonsPM360, minDLon);
        res[i][LON_INDEX] += (indexMin - 1)*360.0;
   }

    int indexPole = -1;
    double avgLon = 0.;
    for (size_t i = 0; i < res.size(); ++i) {
        // detect if node is on/near pole
        if (std::abs(std::abs(res[i][LAT_INDEX]) -  90.) < 1.e-12) {
            indexPole = i;
        }
        else {
            avgLon += res[i][LON_INDEX];
        }
    }
    avgLon /= 3.;

    if (indexPole >= 0) {
        // longitude at the pole is ill defined - we can set it to any
        // value.
        if (avgLon > 180.0) {
            res[indexPole][LON_INDEX] = 270.0;
        }
        else {
            res[indexPole][LON_INDEX] = 90.0;
        }
    }

    return res;
}


std::vector<Vec3> 
Ugrid2D::getEdgePointsRegularized(size_t edgeId) const {

    const size_t* ptIds = this->getEdgePointIds(edgeId);
    const double* p0 = this->getPoint(ptIds[0]);
    const double* p1 = this->getPoint(ptIds[1]);

    std::vector<Vec3> res(2);
    res[0] = Vec3(p0);
    res[1] = Vec3(p1);

    // fix the longitude to minimize the edge length
    double dLon = p1[LON_INDEX] - p0[LON_INDEX];
    double dLonsPM360[] = {std::abs(dLon - 360.), std::abs(dLon), std::abs(dLon + 360.)};
    double* minDLon = std::min_element(&dLonsPM360[0], &dLonsPM360[3]);
    int indexMin = (int) std::distance(dLonsPM360, minDLon);
    res[1][LON_INDEX] += (indexMin - 1)*360.0;

    // fix the latitude to minimize the edge length
    if (std::abs(std::abs(res[0][LAT_INDEX]) - 90.) < 1.e-12) {
        // lon is not well defined at the pole, choose to be the same as the other lon
        res[0][LON_INDEX] = res[1][LON_INDEX];
    }
    if (std::abs(res[1][LAT_INDEX]) == 90.) {
        // lon is not well defined at the pole, choose to be the same as the other lon
        res[1][LON_INDEX] = res[0][LON_INDEX];
    }

    return res;
}

void
Ugrid2D::dumpGridVtk(const std::string& filename) {
    std::ofstream f;
    f.open(filename);
    f << "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    f << "POINTS " << 4 * this->numFaces << " double\n";
    for (size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        const std::vector<Vec3> nodes = this->getFacePointsRegularized(faceId);
        for (const Vec3& node : nodes) {
            f << node << ' ';
        }
        f << '\n';
    }
    f << "CELLS " << this->numFaces << ' ' << 5 * this->numFaces << '\n'; // 2D
    for (size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        // 4 points per face
        f << "4 " << faceId*4 + 0 << ' '
                  << faceId*4 + 1 << ' '
                  << faceId*4 + 2 << ' '
                  << faceId*4 + 3 << '\n';
    }
    f << "CELL_TYPES " << this->numFaces << '\n';
    for (size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        f << "9 ";
        if (faceId % 10 == 0) f << '\n';
    }
    f << '\n';
    f.close();
}


