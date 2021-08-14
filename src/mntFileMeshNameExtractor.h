#include <map>
#include <string>

#ifndef MNT_FILE_MESH_NAME_EXTRACTOR
#define MNT_FILE_MESH_NAME_EXTRACTOR

std::map<std::string, std::string> fileMeshNameExtractor(const std::string& fm) {

    std::map<std::string, std::string> res;

    // extract the filename and the mesh name from "filename:meshname"
    std::size_t columnPosR = fm.rfind(':');
    std::string filename = fm.substr(0, columnPosR);
    std::string meshname = fm.substr(columnPosR + 1, std::string::npos);
    res.insert(std::pair<std::string, std::string>("filename", filename));
    res.insert(std::pair<std::string, std::string>("meshname", meshname));

    return res;
}

std::map<std::string, std::string> fileMeshNameExtractor(const char* fileAndMeshName) {
    std::string fm = std::string(fileAndMeshName);
    return fileMeshNameExtractor(fm);
}

#endif // MNT_FILE_MESH_NAME_EXTRACTOR
