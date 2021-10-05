#include <mntFileMeshNameExtractor.h>

const char FILE_MESH_SEPARATOR = '$';

std::pair<std::string, std::string> fileMeshNameExtractor(const std::string& fm) {
    // extract the filename and the mesh name from "filename:meshname"
    std::size_t columnPosR = fm.rfind(FILE_MESH_SEPARATOR);
    std::string filename = fm.substr(0, columnPosR);
    std::string meshname = fm.substr(columnPosR + 1, std::string::npos);
    return std::pair<std::string, std::string>(filename, meshname);
}

std::pair<std::string, std::string> fileMeshNameExtractor(const char* fileAndMeshName) {
    std::string fm = std::string(fileAndMeshName);
    return fileMeshNameExtractor(fm);
}
