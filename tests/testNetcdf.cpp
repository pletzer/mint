#include <netcdf.h>
#include <string>
#include <iostream>
#undef NDEBUG // turn on asserts
#include <cassert>

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

    ier = nc_close(ncid);
    assert(ier == NC_NOERR);
}

int main(int argc, char** argv) {

	test("${CMAKE_SOURCE_DIR}/data/cs_16.nc", "physics");

    return 0;
}
