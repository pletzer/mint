#include "mntLogger.h"
#include "mntNcFieldRead.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

LIBRARY_API
int mnt_openNc(const char* filename, const char* varname, int* ncid, int* varid) {
  int ier = nc_open(filename, NC_NOWRITE, ncid);
  if (ier != NC_NOERR) {
    return 1;
  }
  ier = nc_inq_varid(*ncid, varname, varid);
  if (ier != NC_NOERR) {
    return 2;
  }
  return 0;
}

LIBRARY_API
int mnt_closeNc(int ncid) {
  int ier = nc_close(ncid);
  return ier;
}

LIBRARY_API
int mnt_ncfieldread_new(NcFieldRead_t** self, int ncid, int varid) {

  int ier;
  std::string msg;

  *self = new NcFieldRead_t();
  (*self)->ncid = ncid;
  (*self)->varid = varid;

  // inquire about the variable
  int ndims = 0;
  ier = nc_inq_varndims(ncid, varid, &ndims);
  if (ier != NC_NOERR) {
    msg = "could not inquire ndims for variable with Id " 
              + std::to_string(varid) + " in netcdf file with Id "
              + std::to_string(ncid) + " ier = " + std::to_string(ier);
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 3;
  }

  int dimIds[NC_MAX_VAR_DIMS];
  ier = nc_inq_vardimid(ncid, varid, dimIds);
  if (ier != NC_NOERR) {
    msg = "could not inquire dim ids for variable with Id " 
              + std::to_string(varid) + " in netcdf file with Id "
              + std::to_string(ncid) + " ier = " + std::to_string(ier);
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 4;
  }

  // get the dimensions and dimension names
  char dimName[NC_MAX_NAME + 1];
  std::size_t dim;
  for (int iDim = 0; iDim < ndims; ++iDim) {
    ier = nc_inq_dim(ncid, dimIds[iDim], dimName, &dim);
    if (ier == NC_NOERR) {
      (*self)->dimNames.push_back(std::string(dimName));
      (*self)->dimSizes.push_back(dim);
    }
  }

  return 0;
}

LIBRARY_API
int mnt_ncfieldread_del(NcFieldRead_t** self) {
  delete *self;
  return 0;
}


LIBRARY_API
int mnt_ncfieldread_getNumDims(NcFieldRead_t** self, int* ndims) {
  *ndims = (int) (*self)->dimSizes.size();
  return 0;
}


LIBRARY_API
int mnt_ncfieldread_getDimName(NcFieldRead_t** self, 
                               int iAxis, char* dimName, int dimNameLen) {
  dimName = strncpy(dimName, (*self)->dimNames[iAxis].c_str(), dimNameLen);
  return 0;
}


LIBRARY_API
int mnt_ncfieldread_getDim(NcFieldRead_t** self, int iAxis, std::size_t* dim) {
  *dim = (*self)->dimSizes[iAxis];
  return 0;
}

LIBRARY_API
int mnt_ncfieldread_data(NcFieldRead_t** self, 
                         double data[]) {
  int ier = nc_get_var_double((*self)->ncid, (*self)->varid, data);
  if (ier != NC_NOERR) {
    std::string msg = "could not read data";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
  }
  return ier;
}

LIBRARY_API
int mnt_ncfieldread_dataSlice(NcFieldRead_t** self, 
                              const std::size_t startInds0[], 
                              const std::size_t counts[], 
                              double data[]) {
  int ier = nc_get_vara_double((*self)->ncid, (*self)->varid, 
                               startInds0, counts, data);
  if (ier != NC_NOERR) {
    std::string msg = "could not read data slice";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
  }
  return ier;
}

