#include <mntLatLon.h>
#include <CmdLineArgParser.h>

int main(int argc, char** argv) {

    CmdLineArgParser args;
    args.setPurpose("Create a lat-lon uniform grid.");
    args.set("-nlat", 2, "Number of latitude cells");
    args.set("-nlon", 4, "Number of longitude cells");
    args.set("-o", std::string("um.nc"), "Output file");

    int nlat = args.get<int>("-nlat");
    int nlon = args.get<int>("-nlon");

    mntLatLon_t *ll;
    mnt_latlon_new(&ll);
    mnt_latlon_del(&ll);
}