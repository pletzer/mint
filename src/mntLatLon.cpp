#include <mntLatLon.h>
#include <netcdf.h>
#include <iostream>


extern "C"
int mnt_latlon_new(LatLon_t** self){
    *self = new LatLon_t();
    (*self)->fillValue = -1.073742e+09;
    (*self)->dLat = 0.;
    (*self)->dLon = 0.;
    return 0;
}

extern "C"
int mnt_latlon_del(LatLon_t** self) {
    delete *self;
    return 0;
}

extern "C"
int mnt_latlon_setNumberOfLatCells(LatLon_t** self, size_t n) {
    (*self)->lats.resize(n + 1);
    return 0;
}

extern "C"
int mnt_latlon_setNumberOfLonCells(LatLon_t** self, size_t n) {
    (*self)->lons.resize(n + 1);
    return 0;
}

extern "C"
int mnt_latlon_build(LatLon_t** self) {
    size_t nLats = (*self)->lats.size();
    size_t nLons = (*self)->lons.size();
    (*self)->dLat = 180.0 / double(nLats - 1);
    (*self)->dLon = 360.0 / double(nLons - 1);
    for (size_t i = 0; i < nLats; ++i) {
        (*self)->lats[i] = -90.0 + i * (*self)->dLat;
    }
    for (size_t i = 0; i < nLons; ++i) {
        (*self)->lons[i] = 0.0 + i * (*self)->dLon;
    }
    return 0;
}

extern "C"
int mnt_latlon_load(LatLon_t** self, const std::string& filename) {

    int ncid, latitude_0_dim, longitude_0_dim;
    size_t numLat0, numLon0;
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
    return 0;
}

extern "C"
int mnt_latlon_dump(LatLon_t** self, const std::string& filename) {

    int ncid, latitude_0_dim, latitude_dim, longitude_0_dim, longitude_dim;
    int latitude_var, longitude_var;
    int dims_lat[1];
    int dims_lon[1];
    int ier = 0;

    size_t numLat = (*self)->lats.size();
    size_t numLon = (*self)->lons.size() - 1; // periodic
    size_t numLat0 = numLat - 1;
    size_t numLon0 = numLon;

    ier = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when creating the file\n";

    // latitude dimension
    ier = nc_def_dim(ncid, "latitude_0", numLat0, &latitude_0_dim);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining dim latitude0\n";
    ier = nc_def_dim(ncid, "latitude", numLat, &latitude_dim);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining dim latitude\n";

    // longitude dimension
    ier = nc_def_dim(ncid, "longitude_0", numLon0, &longitude_0_dim);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining dim longitude0\n";
    ier = nc_def_dim(ncid, "longitude", numLon, &longitude_dim);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining dim longitude\n";

    // variables
    dims_lat[0] = latitude_dim;
    ier = nc_def_var(ncid, "latitude", NC_DOUBLE, 1, dims_lat, &latitude_var);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining var latitude\n";
    ier = nc_put_att_text(ncid, latitude_var, "axis", 1, "Y");
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when adding att axis to var latitude\n";
    ier = nc_put_att_text(ncid, latitude_var, "units", 13, "degrees_north");
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when adding att units to var latitude\n";

    dims_lon[0] = longitude_dim;
    ier = nc_def_var(ncid, "longitude", NC_DOUBLE, 1, dims_lon, &longitude_var);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when defining var longitude\n";
    ier = nc_put_att_text(ncid, longitude_var, "axis", 1, "X");
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when adding att axis to var longitude\n";
    ier = nc_put_att_text(ncid, longitude_var, "units", 12, "degrees_east");
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when adding att units to var longitude\n";

    ier = nc_enddef(ncid);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when ending define mode\n";

    // write data
    ier = nc_put_var_double(ncid, latitude_var, &(*self)->lats[0]);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when writing var latitude\n";
    ier = nc_put_var_double(ncid, longitude_var, &(*self)->lons[0]);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when writing var longitude\n";

    ier = nc_close(ncid);
    if (ier != 0) std::cerr << "mnt_latlon_dump: ERROR when closing file\n";

    return 0;
}


