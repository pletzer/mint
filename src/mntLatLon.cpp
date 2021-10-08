#include "mntLogger.h"
#include <mntLatLon.h>
#include <netcdf.h>
#include <iostream>


LIBRARY_API
int mnt_latlon_new(LatLon_t** self){
    *self = new LatLon_t();
    (*self)->fillValue = -1.073742e+09;
    (*self)->dLat = 0.;
    (*self)->dLon = 0.;
    return 0;
}

LIBRARY_API
int mnt_latlon_del(LatLon_t** self) {
    delete *self;
    return 0;
}

LIBRARY_API
int mnt_latlon_setNumberOfLatCells(LatLon_t** self, std::size_t n) {
    (*self)->lats.resize(n + 1);
    return 0;
}

LIBRARY_API
int mnt_latlon_setNumberOfLonCells(LatLon_t** self, std::size_t n) {
    (*self)->lons.resize(n + 1);
    return 0;
}

LIBRARY_API
int mnt_latlon_build(LatLon_t** self) {
    std::size_t nLats = (*self)->lats.size();
    std::size_t nLons = (*self)->lons.size();
    (*self)->dLat = 180.0 / double(nLats - 1);
    (*self)->dLon = 360.0 / double(nLons - 1);
    for (std::size_t i = 0; i < nLats; ++i) {
        (*self)->lats[i] = -90.0 + i * (*self)->dLat;
    }
    for (std::size_t i = 0; i < nLons; ++i) {
        (*self)->lons[i] = 0.0 + i * (*self)->dLon;
    }
    return 0;
}

LIBRARY_API
int mnt_latlon_load(LatLon_t** self, const std::string& filename) {

    int ncid, latitude_0_dim, longitude_0_dim;
    std::size_t numLat0, numLon0;
    int ier = 0;

    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);

    // latitude dimension
    ier = nc_inq_dimid(ncid, "latitude_0", &latitude_0_dim);
    ier = nc_inq_dimlen(ncid, latitude_0_dim, &numLat0);
    // longitude dimension
    ier = nc_inq_dimid(ncid, "longitude_0", &longitude_0_dim);
    ier = nc_inq_dimlen(ncid, longitude_0_dim, &numLon0);

    ier = nc_close(ncid);

    mnt_latlon_new(self);
    mnt_latlon_setNumberOfLatCells(self, numLat0);
    mnt_latlon_setNumberOfLonCells(self, numLon0);

    mnt_latlon_build(self);
    return ier;
}

LIBRARY_API
int mnt_latlon_dump(LatLon_t** self, const std::string& filename) {

    int ncid, latitude_0_dim, latitude_dim, longitude_0_dim, longitude_dim;
    int latitude_var, longitude_var;
    int dims_lat[1];
    int dims_lon[1];
    int ier = 0;

    std::size_t numLat = (*self)->lats.size();
    std::size_t numLon = (*self)->lons.size() - 1; // periodic
    std::size_t numLat0 = numLat - 1;
    std::size_t numLon0 = numLon;

    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "creating the file");

    // latitude dimension
    ier = nc_def_dim(ncid, "latitude_0", numLat0, &latitude_0_dim);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining dim latitude0");
    ier = nc_def_dim(ncid, "latitude", numLat, &latitude_dim);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining dim latitude");

    // longitude dimension
    ier = nc_def_dim(ncid, "longitude_0", numLon0, &longitude_0_dim);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining dim longitude0");
    ier = nc_def_dim(ncid, "longitude", numLon, &longitude_dim);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining dim longitude");

    // variables
    dims_lat[0] = latitude_dim;
    ier = nc_def_var(ncid, "latitude", NC_DOUBLE, 1, dims_lat, &latitude_var);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining var latitude");
    ier = nc_put_att_text(ncid, latitude_var, "axis", 1, "Y");
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "adding att axis to var latitude");
    ier = nc_put_att_text(ncid, latitude_var, "units", 13, "degrees_north");
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "adding att units to var latitude");

    dims_lon[0] = longitude_dim;
    ier = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, dims_lon, &longitude_var);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "defining var longitude");
    ier = nc_put_att_text(ncid, longitude_var, "axis", 1, "X");
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "adding att axis to var longitude");
    ier = nc_put_att_text(ncid, longitude_var, "units", 12, "degrees_east");
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "adding att units to var longitude");

    ier = nc_enddef(ncid);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "ending define mode");

    // write data
    ier = nc_put_var_double(ncid, latitude_var, &(*self)->lats[0]);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "writing var latitude");
    ier = nc_put_var_double(ncid, longitude_var, &(*self)->lons[0]);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "writing var longitude");

    ier = nc_close(ncid);
    if (ier != NC_NOERR) mntlog::error(__FILE__, __func__, __LINE__, "closing file");

    return 0;
}


