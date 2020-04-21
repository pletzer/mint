#include <netcdf.h>
#include <string>
#include <iostream>
#undef NDEBUG // turn on asserts
#include <cassert>
#include "mntNcDimensions.h"

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

    NcDimensions_t* attrs;
    ier = mnt_ncdimensions_new(&attrs);
    assert(ier == NC_NOERR);

    ier = mnt_ncdimensions_read(&attrs, ncid, varid);
    assert(ier == NC_NOERR);

    int ndims;
    ier = mnt_ncdimensions_getNumDims(&attrs, &ndims);
    assert(ier == NC_NOERR);
    assert(ndims == 3);

    size_t dims[ndims];
    for (int i = 0; i < ndims; ++i) {
        ier = mnt_ncdimensions_get(&attrs, i, &dims[i]);
        assert(ier == NC_NOERR);
    }

    ier = nc_close(ncid);
    assert(ier == NC_NOERR);

    ier = mnt_ncdimensions_print(&attrs);
    assert(ier == NC_NOERR);
}

int main(int argc, char** argv) {

	test("${CMAKE_SOURCE_DIR}/data/lfric_diag_ex.nc", "u2");

    return 0;
}
