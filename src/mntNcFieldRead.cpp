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

  // open the file
  int ier = nc_open(fname.c_str(), NC_NOWRITE, &(*self)->ncid);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not open file " << fname << '\n';
    return 1;
  }

  // get the variable Id
  ier = nc_inq_varid((*self)->ncid, vname.c_str(), &(*self)->varid);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not find variable " << vname << " in file " << fname << '\n';
    return 2;
  }

  // inquire about the variable
  int ndims = 0;
  int dimIds[NC_MAX_VAR_DIMS];
  int natts = 0;
  ier = nc_inq_var((*self)->ncid, (*self)->varid, NULL, NULL, &ndims, dimIds, &natts);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not inquire about variable " << vname << " in file " << fname << '\n';
    return 3;
  }

  // get the dimensions and dimension names
  char dimName[NC_MAX_NAME + 1];
  size_t dim;
  for (int iDim = 0; iDim < ndims; ++iDim) {
    ier = nc_inq_dim((*self)->ncid, dimIds[iDim], dimName, &dim);
    if (ier == NC_NOERR) {
      (*self)->dimNames.push_back(std::string(dimName));
      (*self)->dimSizes.push_back(dim);
    }
  }

  // get the attributes
  char attname[NC_MAX_NAME + 1];
  size_t n;
  nc_type xtype;
  for (int i = 0; i < natts; ++i) {
    ier = nc_inq_attname((*self)->ncid, (*self)->varid, i, attname);
    ier = nc_inq_att((*self)->ncid, (*self)->varid, attname, &xtype, &n);
    if (n == 1 && xtype == NC_DOUBLE) {
      double val;
      ier = nc_get_att_double((*self)->ncid, (*self)->varid, attname, &val);
      (*self)->attDbl.insert(std::pair<std::string, double>(std::string(attname), val));
    }
    else if (n == 1 && xtype == NC_INT) {
      int val;
      ier = nc_get_att_int((*self)->ncid, (*self)->varid, attname, &val);
      (*self)->attInt.insert(std::pair<std::string, int>(std::string(attname), val));
    }
    else if (xtype == NC_CHAR) {
      char val[n + 1];
      ier = nc_get_att_text((*self)->ncid, (*self)->varid, attname, val);
      (*self)->attStr.insert(std::pair<std::string, std::string>(std::string(attname), val));
    }
    else {
      std::cerr << "Warning: unsupported attribute type " << xtype << " of length " << n << '\n';
    }
    if (ier != NC_NOERR) {
      std::cerr << "Warning: failed to read attribute " << attname << " of variable " << vname << '\n';
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
int mnt_ncfieldread_getNumAttsStr(NcFieldRead_t** self,
                                 int* n) {
  *n = (*self)->attStr.size();
  return 0;
}

extern "C"
int mnt_ncfieldread_getNumAttsInt(NcFieldRead_t** self,
                                 int* n) {
  *n = (*self)->attInt.size();
  return 0;
}

extern "C"
int mnt_ncfieldread_getNumAttsDbl(NcFieldRead_t** self,
                                 int* n) {
  *n = (*self)->attDbl.size();
  return 0;
}

extern "C"
int mnt_ncfieldread_getAttsStr(NcFieldRead_t** self,
                              char attNames[], int attNameLen,
                              char attVals[], int attValLen) {

  size_t n = (*self)->attStr.size();
  size_t count = 0;
  for (auto it = (*self)->attStr.cbegin(); it != (*self)->attStr.cend(); ++it) {
    strncpy(&attNames[count*attNameLen], it->first.c_str(), attNameLen);
    strncpy(&attVals[count*attValLen], it->second.c_str(), attValLen);
    count++;
  }
  return 0;
}

extern "C"
int mnt_ncfieldread_getAttsInt(NcFieldRead_t** self,
                              char attNames[], int attNameLen,
                              int attVals[]) {

  size_t n = (*self)->attInt.size();
  size_t count = 0;
  for (auto it = (*self)->attInt.cbegin(); it != (*self)->attInt.cend(); ++it) {
    strncpy(&attNames[count*attNameLen], it->first.c_str(), attNameLen);
    attVals[count] = it->second;
    count++;
  }
  return 0;
}

extern "C"
int mnt_ncfieldread_getAttsDbl(NcFieldRead_t** self,
                              char attNames[], int attNameLen,
                              double attVals[]) {


  size_t n = (*self)->attDbl.size();
  size_t count = 0;
  for (auto it = (*self)->attDbl.cbegin(); it != (*self)->attDbl.cend(); ++it) {
    strncpy(&attNames[count*attNameLen], it->first.c_str(), attNameLen);
    attVals[count] = it->second;
    count++;
  }
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

