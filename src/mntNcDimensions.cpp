#include "mntNcDimensions.h"
#include <netcdf.h>
#include <iostream>

LIBRARY_API
int mnt_ncdimensions_new(NcDimensions_t** self) {
  *self = new NcDimensions_t();
  return 0;
}

LIBRARY_API
int mnt_ncdimensions_del(NcDimensions_t** self) {
  delete *self;
  return 0;
}

LIBRARY_API
int mnt_ncdimensions_read(NcDimensions_t** self, int ncid, int varid) {
  int ier;
  int ndims;

  ier = nc_inq_varndims(ncid, varid, &ndims);
  if (ier != NC_NOERR) return ier;

  (*self)->dims.resize(ndims);
  std::vector<int> dimIds(ndims);
  ier =  nc_inq_vardimid(ncid, varid, &dimIds[0]);
  if (ier != NC_NOERR) return ier;

  for (int i = 0; i < ndims; ++i) {
    ier = nc_inq_dimlen(ncid, dimIds[i], &(*self)->dims[i]);
  }

  return 0;
}

LIBRARY_API
int mnt_ncdimensions_getNumDims(NcDimensions_t** self, int* ndims) {
  *ndims = (int) (*self)->dims.size();
  return 0;
}

LIBRARY_API
int mnt_ncdimensions_get(NcDimensions_t** self, int i, std::size_t* len) {
  *len = (*self)->dims[i];
  return 0;
}

LIBRARY_API
int mnt_ncdimensions_print(NcDimensions_t** self) {
  for (std::size_t i = 0; i < (*self)->dims.size(); ++i) {
    std::cout << i << " size: " << (*self)->dims[i] << '\n';
  }
  return 0;
}

