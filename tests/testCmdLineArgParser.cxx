#include <iostream>
#include <mntCmdLineArgParser.h>


int main(int argc, char** argv) {

    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename$meshname\"");
    args.set("-v", std::string(""), "Specify edge staggered field variable name in source UGRID file, varname[@filename$meshname]");
    args.set("-P", 0.0, "Specify the periodicity length in longitudes (default is non-periodic)");
    args.set("-d", std::string(""), "UGRID destination grid file name");
    args.set("-w", std::string(""), "Write interpolation weights to file");
    args.set("-W", std::string(""), "Load interpolation weights from file");
    args.set("-o", std::string(""), "Specify output VTK file where regridded edge data are saved");
    args.set("-O", std::string(""), "Specify output 2D UGRID file where regridded edge data are saved");
    args.set("-S", 1, "Set to zero to disable source grid regularization, -S 0 is required for uniform lon-lat grid");
    args.set("-D", 1, "Set to zero to disable destination grid regularization, -S 0 is required for uniform lon-lat grid");
    args.set("-N", 128, "Average number of cells per bucket");
    args.set("-debug", 1, "0=no checks, 1=print outside segments, 2=save outside segments");
    args.set("-verbose", false, "Turn on verbosity");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    args.print();

    return 0;
}

