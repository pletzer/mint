#include <CmdLineArgParser.h>

#ifndef MNT_CMD_LINE_ARG_PARSER
#define MNT_CMD_LINE_ARG_PARSER

/**
 * Constructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_new(CmdLineArgParser** self);

/**
 * Destructor
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_del(CmdLineArgParser** self);

/**
 * Set double command line argument
 * @param name name of the command line argument
 * @param def_value default value
 * @param help help message
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_setdouble(CmdLineArgParser** self, const char* name, double def_value, const char* help);

/**
 * Set integer command line argument
 * @param name name of the command line argument
 * @param def_value default value
 * @param help help message
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_setint(CmdLineArgParser** self, const char* name, int def_value, const char* help);

/**
 * Set string command line argument
 * @param name name of the command line argument
 * @param def_value default value
 * @param help help message
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_setstring(CmdLineArgParser** self, const char* name, const char* def_value, const char* help);

/**
 * Parse the command line arguments
 * @param nargs number of arguments
 * @param args_lengths length of each command line argument
 * @param args list of command line arguments
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_parse(CmdLineArgParser** self, int nargs, const int args_lengths[], char** args);

/**
 * Get double command line argument
 * @param name name of the command line argument
 * @param val returned value
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_getdouble(CmdLineArgParser** self, const char* name, double* val);

/**
 * Get integer command line argument
 * @param name name of the command line argument
 * @param val returned value
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_getint(CmdLineArgParser** self, const char* name, int* val);

/**
 * Get string command line argument
 * @param name name of the command line argument
 * @param val string will be coped into val, val must be have been allocated (size >= n + 1)
 * @param n number of characters in the string, excluding '\0' (out)
 * @return error code (0 is OK)
 * @note a termination character '\0' will be added
 */
extern "C"
int mnt_cmdlineargparser_getstring(CmdLineArgParser** self, const char* name, char* val, int* n);



#endif // MNT_CMD_LINE_ARG_PARSER
