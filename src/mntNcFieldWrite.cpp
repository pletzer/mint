#include "mntLogger.h"
#include "mntNcFieldWrite.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

LIBRARY_API
int mnt_ncfieldwrite_new(NcFieldWrite_t** self,
                        const char* fileName, int fileNameLen, 
                        const char* varName, int varNameLen,
                        int append) {
  std::string msg;

  *self = new NcFieldWrite_t();
  (*self)->ncid = -1;
  (*self)->varid = -1;
  (*self)->defined = false;
  (*self)->append = false;
  (*self)->varName = std::string(varName, varNameLen);

  if (append == 1) {
    (*self)->append = true;
  }

  std::string fname = std::string(fileName, fileNameLen);
  std::string vname = std::string(varName, varNameLen);

  if ((*self)->append) {
    // append data to existing file
    int ier = nc_open(fname.c_str(), NC_WRITE, &(*self)->ncid);
    if (ier != NC_NOERR) {
      msg = "could not open file " + fname + " in write access mode";
      mntlog::error(__FILE__, __func__, __LINE__, msg);
      return 1;
    }

    // load the attributes and dimensions if in append mode
    ier = mnt_ncfieldwrite_inquire(self);
    if (ier != NC_NOERR) {
      msg = "failed to inquire the attributes and dimensions for variable " 
            + std::string(varName) + " in file " + fname;
      mntlog::error(__FILE__, __func__, __LINE__, msg);
      return 2;
    }

  }
  else {
    // create a new file
    int ier = nc_create(fname.c_str(), NC_CLOBBER, &(*self)->ncid);
    if (ier != NC_NOERR) {
      msg = "could not create file " + fname;
      mntlog::error(__FILE__, __func__, __LINE__, msg);
      return 1;
    }
  }

  return 0;
}

LIBRARY_API
int mnt_ncfieldwrite_del(NcFieldWrite_t** self) {
  int ier = nc_close((*self)->ncid);
  delete *self;
  return ier;
}


LIBRARY_API
int mnt_ncfieldwrite_setNumDims(NcFieldWrite_t** self, int ndims) {
  (*self)->dimSizes.resize(ndims);
  (*self)->dimNames.resize(ndims);
  return 0;
}


LIBRARY_API
int mnt_ncfieldwrite_setDim(NcFieldWrite_t** self, 
                            int iAxis, 
                            const char* dimName, int dimNameLen,
                            std::size_t dim) {
  (*self)->dimNames[iAxis] = std::string(dimName, dimNameLen);
  (*self)->dimSizes[iAxis] = dim;
  return 0;
}

LIBRARY_API
int mnt_ncfieldwrite_data(NcFieldWrite_t** self, 
                          const double data[]) {

  if (!(*self)->defined) {
    int stat = mnt_ncfieldwrite_define(self);
    if (stat != NC_NOERR) {
       return stat;
    }
  }

  int ier = nc_put_var_double((*self)->ncid, (*self)->varid, data);
  return ier;
}

LIBRARY_API
int mnt_ncfieldwrite_dataSlice(NcFieldWrite_t** self, 
                               const std::size_t startInds0[], 
                               const std::size_t counts[], 
                               const double data[]) {
  if (!(*self)->defined) {
    int stat = mnt_ncfieldwrite_define(self);
    if (stat != NC_NOERR) {
       return stat;
    }
  }

  int ier = nc_put_vara_double((*self)->ncid, (*self)->varid, 
                               startInds0, counts, data);
  return ier;
}

LIBRARY_API
int mnt_ncfieldwrite_define(NcFieldWrite_t** self) {

  std::string msg;

  if ((*self)->defined) {
    // nothing to do
    return 0;
  }

  if ((*self)->append) {
    msg = "can define only in non-append mode\n";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 1;
  }  

  int ndims = (int) (*self)->dimSizes.size();
  std::vector<int> dimIds(ndims);
  int ier;
  for (std::size_t i = 0; i < ndims; ++i) {
    ier = nc_def_dim((*self)->ncid, (*self)->dimNames[i].c_str(),
                                    (*self)->dimSizes[i], &dimIds[i]);
    if (ier != NC_NOERR) {
      msg = "could not define dimension " 
                + (*self)->dimNames[i] + " = " 
                + std::to_string((*self)->dimSizes[i]);
      mntlog::error(__FILE__, __func__, __LINE__, msg);
      return 2;
    }
  }

  // define the variable
  ier = nc_def_var((*self)->ncid, (*self)->varName.c_str(), NC_DOUBLE, ndims, 
                   &dimIds[0], &(*self)->varid);
    if (ier != NC_NOERR) {
      msg = "could not define variable " + (*self)->varName;
      mntlog::error(__FILE__, __func__, __LINE__, msg);
      return 3;
    }

  (*self)->defined = true;
  ier = nc_enddef((*self)->ncid);

  return ier;
}

LIBRARY_API
int mnt_ncfieldwrite_inquire(NcFieldWrite_t** self) {

  int ier;
  std::string msg;

  if (!(*self)->append) {
    msg = "can inquire only in append mode\n";
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 1;
  }

  // get the variable Id
  ier = nc_inq_varid((*self)->ncid, (*self)->varName.c_str(), &(*self)->varid);
  if (ier != NC_NOERR) {
    msg = "could not find variable " + (*self)->varName;
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 2;
  }

  // enquire about the variable
  int ndims = 0;
  int dimIds[NC_MAX_VAR_DIMS];
  int natts = 0;
  ier = nc_inq_var((*self)->ncid, (*self)->varid, NULL, NULL, &ndims, dimIds, &natts);
  if (ier != NC_NOERR) {
    msg = "could not inquire about variable " + (*self)->varName;
    mntlog::error(__FILE__, __func__, __LINE__, msg);
    return 3;
  }

  // get the dimensions and dimension names
  (*self)->dimNames.resize(0);
  (*self)->dimSizes.resize(0);
  char dimName[NC_MAX_NAME + 1];
  std::size_t dim;
  for (int iDim = 0; iDim < ndims; ++iDim) {
    ier = nc_inq_dim((*self)->ncid, dimIds[iDim], dimName, &dim);
    if (ier == NC_NOERR) {
      (*self)->dimNames.push_back(std::string(dimName));
      (*self)->dimSizes.push_back(dim);
    }
  }

  // no need to define the variable in append mode
  (*self)->defined = true;

  return 0;
}
