#include <netcdf.h>
#include <string>
#include <iostream>
#undef NDEBUG // turn on asserts
#include <cassert>
#include "mntNcAttributes.h"

void test(const std::string& filename, const std::string& varname) {

	int ncid, ier;

	std::cout << "Opening file: '" << filename << "'\n";
    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    std::cout << "ncid = " << ncid << '\n';
    assert(ier == NC_NOERR);

    int varid;
	std::cout << "Gettig Id of variable: '" << varname << "'\n";
    ier = nc_inq_varid(ncid, varname.c_str(), &varid);
    std::cout << "varid = " << varid << '\n';
    assert(ier == NC_NOERR);

    NcAttributes_t* attrs;
    ier = mnt_ncattributes_new(&attrs);
    assert(ier == NC_NOERR);

    ier = mnt_ncattributes_read(&attrs, ncid, varid);
    assert(ier == NC_NOERR);

    ier = nc_close(ncid);
    assert(ier == NC_NOERR);

    ier = mnt_ncattributes_print(&attrs);
    assert(ier == NC_NOERR);

    assert(mnt_ncattributes_isIntensive(&attrs) == 0);

    // write the attributes in another file
    int ncid2, varid2;
    ier = nc_create("attributes.nc", NC_CLOBBER|NC_NETCDF4, &ncid2);
    assert(ier == NC_NOERR);

    int* dimIds = NULL;
    ier = nc_def_var(ncid2, varname.c_str(), NC_INT, 0, dimIds, &varid2);
    assert(ier == NC_NOERR);

    ier = mnt_ncattributes_write(&attrs, ncid2, varid2);
    assert(ier == NC_NOERR);

    ier = nc_enddef(ncid2);
    assert(ier == NC_NOERR);

    ier = nc_close(ncid2);
    assert(ier == NC_NOERR);
   

}

int main(int argc, char** argv) {

	test("${CMAKE_SOURCE_DIR}/data/cs_16.nc", "physics");

    return 0;
}
