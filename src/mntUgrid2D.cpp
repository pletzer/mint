#include "mntLogger.h"
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
    bool operator()(const std::pair<std::size_t, std::vector<double> >& x, 
                    const std::pair<std::size_t, std::vector<double> >& y) {
        return (x.second[0] < y.second[0]);
    }
};


std::vector<Vec3> 
Ugrid2D::getFacePoints(std::size_t faceId) const {

    const std::size_t* pointIds = this->getFacePointIds(faceId);
    std::vector<Vec3> res(MNT_NUM_VERTS_PER_QUAD);

    // iterate over the points
    for (std::size_t i = 0; i < MNT_NUM_VERTS_PER_QUAD; ++i) {
        const double* p = this->getPoint(pointIds[i]);
        res[i] = Vec3(p);
    }

    return res;
}


std::vector<Vec3> 
Ugrid2D::getEdgePoints(std::size_t edgeId) const {

    std::vector<Vec3> res;

    // iterate over the 2 points spanning the edge
    for (std::size_t i = 0; i < 2; ++i) {

        // get the point id
        std::size_t pointId = this->edge2Points[i + edgeId*2];

        // get the coordinates of this point
        const double* p = this->getPoint(pointId);

        // add
        res.push_back(Vec3(p));
    }
    return res;
}


int 
Ugrid2D::load(const std::string& filename, const std::string& meshname) {

    int ier = 0;
    int ncid;
    std::string msg;

    // open the file
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    if (ier != NC_NOERR) {
        msg = "cannot open netCDF file \"" + filename + "\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    // mesh variable
    int meshid;
    ier = nc_inq_varid(ncid, meshname.c_str(), &meshid);
    if (ier != NC_NOERR) {
        msg = "cannot find mesh name, netCDF variable \"" + meshname + "\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        nc_close(ncid);
        return 2;
    }

    ier = this->readConnectivityData(ncid, meshid, 
                "face_node_connectivity", this->face2Points);
    if (ier != NC_NOERR) {
        msg = "cannot read the face_node_connectivity of mesh \"" + meshname + "\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        nc_close(ncid);
        return 3;
    }

    ier = this->readConnectivityData(ncid, meshid, 
                "face_edge_connectivity", this->face2Edges);
    /* no longer an error
    if (ier != NC_NOERR) {
        msg = "variable \"" + meshname + "\" does not have attribute \"face_edge_connectivity\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        nc_close(ncid);
        return 4;
    }
    */

    ier = this->readConnectivityData(ncid, meshid, 
                "edge_node_connectivity", this->edge2Points);
    if (ier != NC_NOERR) {
        msg = "cannot read the edge_node_connectivity of mesh \"" + meshname + "\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        nc_close(ncid);
        return 5;
    }

    // read the node coordinates
    ier = this->readPoints(ncid, meshid);
    if (ier != NC_NOERR) {
        msg = "cannot read node coordinates for mesh \"" + meshname + "\"";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        nc_close(ncid);
        return 6;
    }

    std::size_t n = this->face2Points.size();
    this->numFaces = n / MNT_NUM_VERTS_PER_QUAD;
    n = this->edge2Points.size();
    this->numEdges = n / MNT_NUM_VERTS_PER_EDGE;

    return 0;
}

int
Ugrid2D::readConnectivityData(int ncid, int meshid, 
                              const std::string& role,
                              std::vector<std::size_t>& data) {

    int ier;
    std::string msg;

    // get the lengths of the attribute string
    std::size_t len;
    ier = nc_inq_attlen(ncid, meshid, role.c_str(), &len);
    if (ier != NC_NOERR) {
        msg = "cannot inquire attribute length of \"" + role + "\"; " +
        nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 1;
    }

    // read the attribute value, ie the name of the variables we will need to read
    std::string varname(len, ' ');
    ier = nc_get_att_text(ncid, meshid, role.c_str(), &varname[0]);
    if (ier != NC_NOERR) {
        msg = "cannot get attribute text of \"" + role +
              "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 2;
    }

    // fetch the variable Id for this variable name
    int varid;
    ier = nc_inq_varid(ncid, varname.c_str(), &varid);
    if (ier != NC_NOERR) {
        msg = "cannot get variable Id for \"" + varname +
        "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 3;
    }

    // dimensions of the variable to read
    int dimids[2];
    std::size_t n0, n1;

    // read the data
    ier = nc_inq_vardimid(ncid, varid, dimids);
    if (ier != NC_NOERR) {
        msg = "cannot inquire the dimension Ids for \"" + varname +
        "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 4;
    }

    ier = nc_inq_dimlen(ncid, dimids[0], &n0);
    if (ier != NC_NOERR) {
        msg = "cannot inquire dimension 0 for \"" + varname +
        "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 5;
    }

    ier = nc_inq_dimlen(ncid, dimids[1], &n1);
    if (ier != NC_NOERR) {
        msg = "cannot inquire dimension 1 for \"" + varname +
        "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 6;
    }

    int startIndex = 0;
    ier = nc_get_att_int(ncid, varid, "start_index", &startIndex);
    if (ier != NC_NOERR) {
        msg = "cannot get attribute value \"start_index\" of variable \"" +
        varname + "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 8;
    }

    std::vector<unsigned long long> buffer(n0 * n1);
    ier = nc_get_var_ulonglong(ncid, varid, &buffer[0]);
    if (ier != NC_NOERR) {
        msg = "cannot read the values of \"" + varname +
        "\"; " + nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 7;
    }

    // subtract start_index
    data.resize(n0 * n1);
    for (std::size_t i = 0; i < n0 * n1; ++i) {
        data[i] = buffer[i] - startIndex;
    }

    return 0;
}

int
Ugrid2D::readPoints(int ncid, int meshid) {

    int ier;
    std::string msg;

    // get the lengths of the attribute string
    std::size_t len;
    ier = nc_inq_attlen(ncid, meshid, "node_coordinates", &len);
    if (ier != NC_NOERR) {
        msg = "cannot inquire attribute \"node_coordinates\"; " +
        std::string(nc_strerror(ier));
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 10;
    }

    // read the attribute, lists the name of the lon and lat coordinates
    std::string val(len, ' ');
    ier = nc_get_att_text(ncid, meshid, "node_coordinates", &val[0]);
    if (ier != NC_NOERR) {
        msg = "cannot get attribute \"node_coordinates\"; " +
        std::string(nc_strerror(ier));
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 11;
    }

    // val is "varx vary" where var{x,y} are the variable names
    std::size_t n = val.size();
    std::size_t spaceL = val.find(' ');
    std::size_t spaceR = val.rfind(' ');
    if (spaceL >= n) {
        // could not find space
        msg = "node_coordinates attribute \"" + val + "\" should contain a space-separated list of coordinate names";
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 19;
    }
    std::string varx = val.substr(0, spaceL);
    std::string vary = val.substr(spaceR + 1, n - spaceR - 1);

    int varxid, varyid;
    ier = nc_inq_varid(ncid, varx.c_str(), &varxid);
    if (ier != NC_NOERR) {
        msg = "cannot inquire \"" + varx + "\"; " +
        nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 12;
    }
    ier = nc_inq_varid(ncid, vary.c_str(), &varyid);
    if (ier != NC_NOERR) {
        msg = "cannot inquire \"" + vary + "\"; " +
        nc_strerror(ier);
        mntlog::error(__FILE__, __func__, __LINE__, msg);
        return 13;
    }

    int dimids[1];

    this->numPoints = 0;
    this->points.resize(0);
    int varids[] = {varxid, varyid};
    for (int ivar = 0; ivar < 2; ++ivar) {

        int varid = varids[ivar];

        // get the attribute length
        ier = nc_inq_attlen(ncid, varid, "standard_name", &len);
        if (ier != NC_NOERR) {
            msg = "cannot inquire attribute \"standard_name\" of variable with varid = " +
            std::to_string(varid) + "; " + nc_strerror(ier);
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 14;
        }
        std::string var_stdn(len, ' ');

        // read the attribute
        ier = nc_get_att_text(ncid, varid, "standard_name", &var_stdn[0]);
        if (ier != NC_NOERR) {
            msg = "cannot get value of attribute \"standard_name\" of variable with varid = " +
            std::to_string(varid) + "; " + nc_strerror(ier);
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 15;
        }

        // get the dimension
        ier = nc_inq_vardimid(ncid, varid, dimids);
        if (ier != NC_NOERR) {
            msg = "cannot inquire the dimension Ids of variable with varid = " +
            std::to_string(varid) + "; " + nc_strerror(ier);
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 16;
        }

        ier = nc_inq_dimlen(ncid, dimids[0], &this->numPoints);
        if (ier != NC_NOERR) {
            msg = "cannot inquire dimension 0 of variable with varid = " +
            std::to_string(varid) + "; " + nc_strerror(ier);
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 17;
        }

        // allocate/resize
        std::vector<double> data(this->numPoints);
        if (this->points.size() == 0) {
            this->points.resize(this->numPoints * NUM_SPACE_DIMS, 0.0);
        }

        // read the data
        ier = nc_get_var_double(ncid, varid, &data[0]);
        if (ier != NC_NOERR) {
            msg = "could not read variable with varid = " + std::to_string(varid);
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 18;
        }
        
        // associate data to coordinate variable
        std::size_t j;
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
            msg = "unknown coordinate with standard_name \"" + var_stdn;
            mntlog::error(__FILE__, __func__, __LINE__, msg);
            return 19; 
        }
        for (std::size_t i = 0; i < this->numPoints; ++i) {
            this->points[j + NUM_SPACE_DIMS*i] = data[i];
        }
    }

    return 0;
}


std::vector<Vec3> 
Ugrid2D::getEdgePointsRegularized(std::size_t edgeId) const {

    const std::size_t* ptIds = this->getEdgePointIds(edgeId);
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
    f << "POINTS " << MNT_NUM_VERTS_PER_QUAD * this->numFaces << " double\n";
    for (std::size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        const std::vector<Vec3> nodes = this->getFacePoints(faceId);
        for (const Vec3& node : nodes) {
            f << node << ' ';
        }
        f << '\n';
    }
    f << "CELLS " << this->numFaces << ' ' << 5 * this->numFaces << '\n'; // 2D
    for (std::size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        f << MNT_NUM_VERTS_PER_QUAD << faceId*MNT_NUM_VERTS_PER_QUAD + 0 << ' '
                                    << faceId*MNT_NUM_VERTS_PER_QUAD + 1 << ' '
                                    << faceId*MNT_NUM_VERTS_PER_QUAD + 2 << ' '
                                    << faceId*MNT_NUM_VERTS_PER_QUAD + 3 << '\n';
    }
    f << "CELL_TYPES " << this->numFaces << '\n';
    for (std::size_t faceId = 0; faceId < this->numFaces; ++faceId) {
        f << "9 ";
        if (faceId % 10 == 0) f << '\n';
    }
    f << '\n';
    f.close();
}


