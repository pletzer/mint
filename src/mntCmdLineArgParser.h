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
 * Set boolean command line argument
 * @param name name of the command line argument
 * @param def_value default value, 0 = false
 * @param help help message
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_setbool(CmdLineArgParser** self, const char* name, int def_value, const char* help);

/**
 * Parse the command line arguments
 * @param nargs number of arguments
 * @param n max length of element of args
 * @param args list of command line arguments as a contiguous array
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_parse(CmdLineArgParser** self, int nargs, int n, char* args);

/**
 * Print help message
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_help(CmdLineArgParser** self);


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
 * @param n length of alocated val
 * @param val string will be coped into val, val must be have been allocated (size >= n + 1)
 * @return error code (0 is OK)
 * @note termination character '\0' will be added to remaining characters that were allocated to val
 */
extern "C"
int mnt_cmdlineargparser_getstring(CmdLineArgParser** self, const char* name, int n, char* val);

/**
 * Get boolean command line argument
 * @param name name of the command line argument
 * @param val returned value, 0 = false
 * @return error code (0 is OK)
 */
extern "C"
int mnt_cmdlineargparser_getbool(CmdLineArgParser** self, const char* name, int* val);


#endif // MNT_CMD_LINE_ARG_PARSER
