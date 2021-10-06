#include "mntLIBRARY_API.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#ifndef MNT_LOGGER
#define MNT_LOGGER

static std::vector<std::string> MNT_LOGS;

namespace mntlog {

    void logging(const std::string& severity, const char* file, const char* function, int lineno, const char* msg) {
        std::string message = severity + ' ' + std::string(file) + 
                              std::string(" in function ") + std::string(function) + 
                              std::string(" (line ") + std::to_string(lineno) + 
                              std::string("): ") + std::string(msg) + '\n';
        MNT_LOGS.push_back(message);
    }

    /** 
     * Log info message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     * */
    void info(const char* file, const char* function, int lineno, const char* msg) {
        logging("info    === ", file, function, lineno, msg);
    }

    /** 
     * Log warning message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     * */
    void warn(const char* file, const char* function, int lineno, const char* msg) {
        logging("Warning === ", file, function, lineno, msg);
    }

    /** 
     * Log error message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     * */
    void error(const char* file, const char* function, int lineno, const char* msg) {
        logging("ERROR   === ", file, function, lineno, msg);
    }

}

/** 
 * Print log messages
 */
LIBRARY_API
void mntPrintLogMessages() {
    for (auto log = MNT_LOGS.begin(); log != MNT_LOGS.end(); ++log) {
        std::cout << *log;
    }
}

/**
 * Write log messages to file (callable from Fortran)
 * @param filename file name
 * @param filename_len number of characters in above filename, excluding '\0'
 * @note this function is callable from Fortran
 */
LIBRARY_API
void mntWriteLogMessages(const char* filename, int filename_len) {
    std::string fn = std::string(filename, filename_len);
    std::ofstream f;
    f.open(filename);
    for (auto log = MNT_LOGS.begin(); log != MNT_LOGS.end(); ++log) {
        f << *log;
    }
    f.close();
}


#endif // MNT_LOGGER