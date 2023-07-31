#include <mntLogger.h>
#include <ctime>


static char STRTIME[32];

void mntlog::logging(const std::string& severity, 
                     const char* file, 
                     const char* function, 
                     int lineno, 
                     const std::string& msg) {

    std::time_t now = std::time(nullptr);

#ifdef __STDC_LIB_EXT1__
    struct tm buf;
    std::strftime(STRTIME, sizeof(STRTIME), "[%c ", localtime_r(&now, &buf));
#else
    // unsafe version
    std::strftime(STRTIME, sizeof(STRTIME), "[%c ", std::localtime(&now));
#endif

    std::string message = std::string(STRTIME) + severity + ' ' + std::string(file) + 
                          std::string(" in function ") + std::string(function) + 
                          std::string(" (line ") + std::to_string(lineno) + 
                          std::string("): ") + std::string(msg) + '\n';
    MNT_LOGS.push_back(message);
}

void mntlog::info(const char* file, const char* function, int lineno,
                  const std::string& msg) {
    logging("info    ] ", file, function, lineno, msg);
}

void mntlog::warn(const char* file, const char* function, int lineno,
                  const std::string& msg) {
    logging("Warning ] ", file, function, lineno, msg);
}

void mntlog::error(const char* file, const char* function, int lineno,
                   const std::string& msg) {
    logging("ERROR   ] ", file, function, lineno, msg);
}

LIBRARY_API
void mnt_printLogMessages() {
    for (auto log = MNT_LOGS.begin(); log != MNT_LOGS.end(); ++log) {
        std::cout << *log;
    }
}

LIBRARY_API
void mnt_writeLogMessages(const char* filename, std::size_t filename_len) {
    std::string fn = std::string(filename, filename_len);
    std::ofstream f;
    f.open(filename);
    for (auto log = MNT_LOGS.begin(); log != MNT_LOGS.end(); ++log) {
        f << *log;
    }
    f.close();
}
