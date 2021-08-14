#include <map>
#include <string>

#ifndef MNT_FILE_MESH_NAME_EXTRACTOR
#define MNT_FILE_MESH_NAME_EXTRACTOR

/**
 * Extract the file and mesh name from string filename:meshname
 * @param fm input string
 * @return map with entries "filename" and "meshname"
 * @note there is no error checking, the retured filename and/or meshname may be empty
 */
std::map<std::string, std::string> fileMeshNameExtractor(const std::string& fm);
std::map<std::string, std::string> fileMeshNameExtractor(const char* fileAndMeshName);

#endif // MNT_FILE_MESH_NAME_EXTRACTOR
