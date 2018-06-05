#include <mntLatLon.h>
#include <CmdLineArgParser.h>
#include <iostream>

int main(int argc, char** argv) {

    CmdLineArgParser args;
    args.setPurpose("Create a lat-lon uniform grid.");
    args.set("-nlat", 2, "Number of latitude cells");
    args.set("-nlon", 4, "Number of longitude cells");
    args.set("-o", std::string("um.nc"), "Output file");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        int nlat = args.get<int>("-nlat");
        int nlon = args.get<int>("-nlon");

        LatLon_t *lalo;
        mnt_latlon_new(&lalo);
        mnt_latlon_setNumberOfLatCells(&lalo, (size_t) args.get<int>("-nlat"));
        mnt_latlon_setNumberOfLonCells(&lalo, (size_t) args.get<int>("-nlon"));
        mnt_latlon_build(&lalo);
        mnt_latlon_dump(&lalo, args.get<std::string>("-o"));
        mnt_latlon_del(&lalo);
    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}