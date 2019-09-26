#include <mntUgrid2D.h>
#include <CmdLineArgParser.h>
#include <cstring>
#undef NDEBUG // turn on asserts


int getFileMeshName(const std::string& fileMeshName, std::string& fileName, std::string& meshName) {
    fileName = "";
    meshName = "";
    size_t columnPosL = fileMeshName.find(':');
    size_t columnPosR = fileMeshName.rfind(':');
    if (columnPosL == std::string::npos) {
       return 1;
    }
    fileName = fileMeshName.substr(0, columnPosL);
    meshName = fileMeshName.substr(columnPosR + 1, std::string::npos);
    return 0;
}

int main(int argc, char** argv) {

    int ier;
    CmdLineArgParser args;
    args.setPurpose("Read a UGRID2D NetCDF file.");
    args.set("-s", std::string(""), "UGRID source grid file and mesh name, specified as \"filename:meshname\"");
    args.set("-vtk", std::string(""), "Output grid VTK file");


    bool success = args.parse(argc, argv);
    bool help = args.get<bool>("-h");

    if (success && !help) {
        std::string srcFile = args.get<std::string>("-s");
        if (srcFile.size() == 0) {
            std::cerr << "ERROR: must specify a source grid file (-s)\n";
            return 1;            
        }
        std::string fileName, meshName;
        int stat = getFileMeshName(srcFile, fileName, meshName);
        if (stat > 0) {
            std::cerr << "ERROR: source grid file and. mesh name must be in the format \"filename:meshname\"\n";
            return 2;
        }

        Ugrid2D ugr;
        ier = ugr.load(fileName, meshName);
        assert(ier == 0);

        std::string vtkOutput = args.get<std::string>("-vtk");
        if (vtkOutput.size() != 0) {
            ugr.dumpGridVtk(vtkOutput);
        }

    }

    return 0;
}   
