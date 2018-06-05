#include <mntRegridEdges.h>
#include <mntGrid.h>
#include <CmdLineArgParser.h>
#include <vtkUnstructuredGrid.h>
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

        // read/build the src grid
        mntGrid_t* srcGrid = NULL;
        mnt_grid_new(&srcGrid);
        mnt_grid_load(&srcGrid, srcFile.c_str());

        // read/build the dst grid
        mntGrid_t* dstGrid = NULL;
        mnt_grid_new(&dstGrid);
        mnt_grid_load(&dstGrid, dstFile.c_str());

        vtkUnstructuredGrid* sg;
        vtkUnstructuredGrid* dg;

        // compute the interpolation weights
        mntRegridEdges_t* rge = NULL;
        mnt_regridedges_new(&rge);
        mnt_regridedges_setSrcGrid(&rge, sg);
        mnt_regridedges_setDstGrid(&rge, dg);
        mnt_regridedges_build(&rge);

        // cleanup
        mnt_grid_del(&dstGrid);
        mnt_grid_del(&srcGrid);
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