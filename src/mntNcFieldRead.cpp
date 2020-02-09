#include "mntNcFieldRead.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

extern "C"
int mnt_ncfieldread_new(NcFieldRead_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen) {

  *self = new NcFieldRead_t();
  (*self)->ncid = -1;
  (*self)->varid = -1;

  std::string fname = std::string(fileName, fileNameLen);
  std::string vname = std::string(varName, varNameLen);

  int ier = nc_open(fname.c_str(), NC_NOWRITE, &(*self)->ncid);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not open file " << fname << '\n';
    return 1;
  }

  ier = nc_inq_varid((*self)->ncid, vname.c_str(), &(*self)->varid);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not find variable " << vname << " in file " << fname << '\n';
    return 2;
  }

  int ndims;
  int dimIds[NC_MAX_VAR_DIMS];
  int natts;
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not inquire about variable " << vname << " in file " << fname << '\n';
    return 3;
  }

  char dimName[NC_MAX_NAME + 1];
  size_t n;
  for (int iDim = 0; iDim < ndims; ++iDim) {
    ier = nc_inq_dim((*self)->ncid, dimIds[iDim], dimName, &n);
    if (ier == NC_NOERR) {
      (*self)->dimNames.push_back(std::string(dimName));
      (*self)->dimSizes.push_back(n);
    }
  }

  return 0;
}

extern "C"
int mnt_ncfieldread_del(NcFieldRead_t** self) {
  int ier = nc_close((*self)->ncid);
  delete *self;
  return ier;
}


extern "C"
int mnt_ncfieldread_getNumDims(NcFieldRead_t** self, int* ndims) {
  *ndims = (*self)->dimSizes.size();
  return 0;
}


extern "C"
int mnt_ncfieldread_getDimName(NcFieldRead_t** self, int iAxis, char* dimName, int dimNameLen) {
  dimName = strncpy(dimName, (*self)->dimNames[iAxis].c_str(), dimNameLen);
  return 0;
}


extern "C"
int mnt_ncfieldread_getDim(NcFieldRead_t** self, int iAxis, size_t* dim) {
  *dim = (*self)->dimSizes[iAxis];
  return 0;
}

extern "C"
int mnt_ncfieldread_readData(NcFieldRead_t** self, 
                             double data[]) {

  int ier = nc_get_var_double((*self)->ncid, (*self)->varid, data);
  return ier;
}

extern "C"
int mnt_ncfieldread_readDataSlice(NcFieldRead_t** self, 
                                  const size_t* startInds0, 
                                  const size_t* counts, 
                                  double data[]) {

  int ier = nc_get_vara_double((*self)->ncid, (*self)->varid, startInds0, counts, data);
  return ier;
}

