#include <mntRegridEdges.h>
#include <CmdLineArgParser.h>
#include <iostream>

int main(int argc, char** argv) {

    CmdLineArgParser args;
    args.setPurpose("Regrid an edge centred field.");
    args.set("-s", std::string(""), "Source grid file in VTK format");
    args.set("-d", std::string(""), "Destination grid file in VTK format");

    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-s");
        std::string dstFile = args.get<std::string>("-d");

        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-s)\n";
            return 1;
        }
        if (dstFile.size() == 0) {
            std::cerr << "ERROR: must specify a destination grid file (-d)\n";
            return 2;
        }

        mntRegridEdges_t* rge = NULL;
        mnt_regridedges_new(&rge);

        // cleanup
        mnt_regridedges_del(&rge);

    }
    else if (help) {
        args.help();
    }
    else {
        std::cerr << "ERROR when parsing command line arguments\n";
    }

    return 0;
}