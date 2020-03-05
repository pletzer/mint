#include "mntNcFieldWrite.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

extern "C"
int mnt_ncfieldwrite_new(NcFieldwrite_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen) {

  *self = new NcFieldwrite_t();
  (*self)->ncid = -1;
  (*self)->varid = -1;
  (*self)->defined = false;
  (*self)->varName = std::string(varName, varNameLen);

  std::string fname = std::string(fileName, fileNameLen);
  std::string vname = std::string(varName, varNameLen);

  // open the file
  int ier = nc_create(fname.c_str(), NC_CLOBBER, &(*self)->ncid);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not create file " << fname << '\n';
    return 1;
  }

  return 0;
}

extern "C"
int mnt_ncfieldwrite_del(NcFieldwrite_t** self) {
  int ier = nc_close((*self)->ncid);
  delete *self;
  return ier;
}


extern "C"
int mnt_ncfieldwrite_setNumDims(NcFieldwrite_t** self, int ndims) {
  (*self)->dimSizes.resize(ndims);
  (*self)->dimNames.resize(ndims);
  return 0;
}


extern "C"
int mnt_ncfieldwrite_setDim(NcFieldwrite_t** self, 
                            int iAxis, 
                            const char* dimName, int dimNameLen,
                            size_t dim) {
  (*self)->dimNames[iAxis] = std::string(dimName, dimNameLen);
  (*self)->dimSizes[iAxis] = dim;
  return 0;
}

extern "C"
int mnt_ncfieldwrite_setAttStr(NcFieldwrite_t** self, 
                               const char* attName, int attNameLen,
                               const char* attVal, int attValLen) {
  std::string nm = std::string(attName, attNameLen);
  std::string val = std::string(attVal, attValLen);
  (*self)->attStr[nm] = val;

  return 0;
}

extern "C"
int mnt_ncfieldwrite_setAttInt(NcFieldwrite_t** self,
                              const char* attName, int attNameLen,
                              int attVal) {
  std::string nm = std::string(attName, attNameLen);
  (*self)->attInt[nm] = attVal;
  return 0;
}

extern "C"
int mnt_ncfieldwrite_setAttDbl(NcFieldwrite_t** self,
                              const char* attName, int attNameLen,
                              double attVal) {
  std::string nm = std::string(attName, attNameLen);
  (*self)->attDbl[nm] = attVal;
  return 0;
}

extern "C"
int mnt_ncfieldwrite_data(NcFieldwrite_t** self, 
                          const double data[]) {

  if (!(*self)->defined) {
    int stat = mnt_ncfieldwrite_define(self);
  }
  int ier = nc_put_var_double((*self)->ncid, (*self)->varid, data);
  return ier;
}

extern "C"
int mnt_ncfieldwrite_dataSlice(NcFieldwrite_t** self, 
                               const size_t startInds0[], 
                               const size_t counts[], 
                               const double data[]) {
  if (!(*self)->defined) {
    int stat = mnt_ncfieldwrite_define(self);
  }
  int ier = nc_put_vara_double((*self)->ncid, (*self)->varid, 
                               startInds0, counts, data);
  return ier;
}

extern "C"
int mnt_ncfieldwrite_define(NcFieldwrite_t** self) {

  size_t ndims = (*self)->dimSizes.size();
  std::vector<int> dimIds(ndims);
  int ier;
  for (int i = 0; i < ndims; ++i) {
    ier = nc_def_dim((*self)->ncid, (*self)->dimNames[i].c_str(),
                                    (*self)->dimSizes[i], &dimIds[i]);
  }

  // define the variable
  ier = nc_def_var((*self)->ncid, (*self)->varName.c_str(), NC_DOUBLE, ndims, 
                   &dimIds[0], &(*self)->varid);

  // add the attributes
  for (auto it = (*self)->attStr.begin(); it != (*self)->attStr.end(); ++it) {
    ier = nc_put_att_text((*self)->ncid, (*self)->varid, it->first.c_str(), it->second.size(), it->second.c_str());
  }
  for (auto it = (*self)->attInt.begin(); it != (*self)->attInt.end(); ++it) {
    ier = nc_put_att_int((*self)->ncid, (*self)->varid, it->first.c_str(), NC_INT, 1, &it->second);
  }
  for (auto it = (*self)->attDbl.begin(); it != (*self)->attDbl.end(); ++it) {
    ier = nc_put_att_double((*self)->ncid, (*self)->varid, it->first.c_str(), NC_DOUBLE, 1, &it->second);
  }

}

