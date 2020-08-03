#include "mntNcAttributes.h"
#include <cstdio>
#include <cstring>
#include <netcdf.h>
#include <iostream>

extern "C"
int mnt_ncattributes_new(NcAttributes_t** self) {
  *self = new NcAttributes_t();
  return 0;
}

extern "C"
int mnt_ncattributes_del(NcAttributes_t** self) {
  (*self)->attStr.clear();
  (*self)->attInt.clear();
  (*self)->attDbl.clear();
  delete *self;
  return 0;
}

extern "C"
int mnt_ncattributes_read(NcAttributes_t** self, int ncid, int varid) {

  int ier;

  // inquire about the variable
  int natts = 0;
  ier = nc_inq_var(ncid, varid, NULL, NULL, NULL, NULL, &natts);
  if (ier != NC_NOERR) {
    std::cerr << "ERROR: could not inquire about variable with Id " << varid << '\n';
    return 3;
  }

  // get the attributes
  char attname[NC_MAX_NAME + 1];
  size_t n;
  nc_type xtype;
  for (int i = 0; i < natts; ++i) {
    ier = nc_inq_attname(ncid, varid, i, attname);
    ier = nc_inq_att(ncid, varid, attname, &xtype, &n);
    if (n == 1 && xtype == NC_DOUBLE) {
      double val;
      ier = nc_get_att_double(ncid, varid, attname, &val);
      (*self)->attDbl.insert(std::pair<std::string, double>(std::string(attname), val));
    }
    else if (n == 1 && xtype == NC_INT) {
      int val;
      ier = nc_get_att_int(ncid, varid, attname, &val);
      (*self)->attInt.insert(std::pair<std::string, int>(std::string(attname), val));
    }
    else if (xtype == NC_CHAR) {
      char val[n + 1];
      ier = nc_get_att_text(ncid, varid, attname, val);
      (*self)->attStr.insert(std::pair<std::string, std::string>(std::string(attname), std::string(val, 0, n)));
    }
    else {
      std::cerr << "Warning: unsupported attribute type " << xtype << " of length " << n << '\n';
    }
    if (ier != NC_NOERR) {
      std::cerr << "Warning: failed to read attribute " << attname << " of variable with Id " << varid << '\n';
    }
  }

  return 0;
}

extern "C"
int mnt_ncattributes_write(NcAttributes_t** self, int ncid, int varid) {

  int ier;

  for (auto it = (*self)->attStr.cbegin(); it != (*self)->attStr.cend(); ++it) {
    ier = nc_put_att_text(ncid, varid, 
                          it->first.c_str(), it->second.size(), it->second.c_str());
    if (ier != NC_NOERR) {
      std::cerr << "ERROR: could not put attribute " 
                << it->first << " = " << it->second << '\n';
      return 4;
    }

  }
  for (auto it = (*self)->attInt.cbegin(); it != (*self)->attInt.cend(); ++it) {
    ier = nc_put_att_int(ncid, varid, 
                         it->first.c_str(), NC_INT, 1, &it->second);
    if (ier != NC_NOERR) {
      std::cerr << "ERROR: could not put attribute " 
                << it->first << " = " << it->second << '\n';
      return 5;
    }
  }
  for (auto it = (*self)->attDbl.cbegin(); it != (*self)->attDbl.cend(); ++it) {
    ier = nc_put_att_double(ncid, varid, 
                            it->first.c_str(), NC_DOUBLE, 1, &it->second);
    if (ier != NC_NOERR) {
      std::cerr << "ERROR: could not put attribute " 
                << it->first << " = " << it->second << '\n';
      return 6;
    }
  }

  return 0;
}

extern "C"
int mnt_ncattributes_isIntensive(NcAttributes_t** self) {
  if ((*self)->attStr["field_methods"] == "intensive") {
    return 1;
  }
  return 0;
}

extern "C"
int mnt_ncattributes_print(NcAttributes_t** self) {

  std::cout << "string attributes:\n";
  for (auto it = (*self)->attStr.cbegin(); it != (*self)->attStr.cend(); ++it) {
    std::cerr << it->first << " -> \"" << it->second << "\"\n";
  }
  std::cout << "int attributes   :\n";
  for (auto it = (*self)->attInt.cbegin(); it != (*self)->attInt.cend(); ++it) {
    std::cerr << it->first << " -> \"" << it->second << "\"\n";
  }
  std::cout << "double attributes:\n";
  for (auto it = (*self)->attDbl.cbegin(); it != (*self)->attDbl.cend(); ++it) {
    std::cerr << it->first << " -> \"" << it->second << "\"\n";
  }

  return 0;
}

