#include "mntLIBRARY_API.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#ifndef MNT_LOGGER
#define MNT_LOGGER

static std::vector<std::string> MNT_LOGS;

namespace mntlog {

    void logging(const std::string& severity, 
    	         const char* file, 
    	         const char* function, 
    	         int lineno, 
    	         const char* msg);

    /** 
     * Log info message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     */
    void info(const char* file, const char* function, int lineno, const char* msg);

    /** 
     * Log warning message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     */
    void warn(const char* file, const char* function, int lineno, const char* msg);
    
    /** 
     * Log error message
     * @param file file name
     * @param function function name
     * @param lineno line number
     * @param msg message
     */
    void error(const char* file, const char* function, int lineno, const char* msg);

}

/** 
 * Print log messages
 */
LIBRARY_API void mnt_printLogMessages();

/**
 * Write log messages to file (callable from Fortran)
 * @param filename file name
 * @param filename_len number of characters in above filename, excluding '\0'
 * @note this function is callable from Fortran
 */
LIBRARY_API void mnt_writeLogMessages(const char* filename, std::size_t filename_len);

#endif // MNT_LOGGER