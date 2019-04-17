#include <netcdf.h>
#include <string>
#undef NDEBUG // turn on asserts
#include <cassert>

void test(const std::string& filename, const std::string& meshname) {

	int ncid, ier;

    ier = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
    assert(ier == NC_NOERR);

    // mesh variable
    int meshid;
    ier = nc_inq_varid(ncid, meshname.c_str(), &meshid);
    assert(ier == NC_NOERR);

    ier = nc_close(ncid);
    assert(ier == NC_NOERR);
}

int main(int argc, char** argv) {

	test("${CMAKE_SOURCE_DIR}/data/cs_16.nc", "physics");

    return 0;
}
