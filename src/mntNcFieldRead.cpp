#include "mntNcFieldRead.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

extern "C"
int mnt_ncfieldread_new(NcFieldRead_t** self, int ncid, int varid) {

  int ier;

  *self = new NcFieldRead_t();
  (*self)->ncid = ncid;
  (*self)->varid = varid;

  // inquire about the variable
  int ndims = 0;
  ier = nc_inq_varndims(ncid, varid, &ndims);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not inquire ndims for variable with Id " 
              << varid << " in netcdf file with Id " << ncid << " ier = " << ier << '\n';
    return 3;
  }

  int dimIds[NC_MAX_VAR_DIMS];
  ier = nc_inq_vardimid(ncid, varid, dimIds);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not inquire dim ids for variable with Id " 
              << varid << " in netcdf file with Id " << ncid << " ier = " << ier << '\n';
    return 4;
  }

  // get the dimensions and dimension names
  char dimName[NC_MAX_NAME + 1];
  size_t dim;
  for (int iDim = 0; iDim < ndims; ++iDim) {
    ier = nc_inq_dim(ncid, dimIds[iDim], dimName, &dim);
    if (ier == NC_NOERR) {
      (*self)->dimNames.push_back(std::string(dimName));
      (*self)->dimSizes.push_back(dim);
    }
  }

  return 0;
}

extern "C"
int mnt_ncfieldread_del(NcFieldRead_t** self) {
  delete *self;
  return 0;
}


extern "C"
int mnt_ncfieldread_getNumDims(NcFieldRead_t** self, int* ndims) {
  *ndims = (int) (*self)->dimSizes.size();
  return 0;
}


extern "C"
int mnt_ncfieldread_getDimName(NcFieldRead_t** self, 
                               int iAxis, char* dimName, int dimNameLen) {
  dimName = strncpy(dimName, (*self)->dimNames[iAxis].c_str(), dimNameLen);
  return 0;
}


extern "C"
int mnt_ncfieldread_getDim(NcFieldRead_t** self, int iAxis, size_t* dim) {
  *dim = (*self)->dimSizes[iAxis];
  return 0;
}

extern "C"
int mnt_ncfieldread_data(NcFieldRead_t** self, 
                         double data[]) {

  int ier = nc_get_var_double((*self)->ncid, (*self)->varid, data);
  return ier;
}

extern "C"
int mnt_ncfieldread_dataSlice(NcFieldRead_t** self, 
                              const size_t startInds0[], 
                              const size_t counts[], 
                              double data[]) {
  int ier = nc_get_vara_double((*self)->ncid, (*self)->varid, 
                               startInds0, counts, data);
  return ier;
}

