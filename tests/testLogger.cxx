#include <mntLogger.h>

int main() {

    mntlog::info(__FILE__, __func__, __LINE__, "an info message");
    mntlog::warn(__FILE__, __func__, __LINE__, "a warning message");
    mntlog::error(__FILE__, __func__, __LINE__, "an error message");

    mntPrintLogMessages();
    std::string filename = "mint_log.txt";
    mntWriteLogMessages(filename.c_str(), (int) filename.size());

    return 0;
}