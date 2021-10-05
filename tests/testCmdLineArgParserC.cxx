#include <mntCmdLineArgParser.h>
#include <cassert>
#include <cstring>
#include <iostream>

int main(int argc, char** argv) {

    CmdLineArgParser* parser;
    int ier;

    // create 

    ier = mnt_cmdlineargparser_new(&parser);
    assert(ier == 0);

    // set the command line options

    ier = mnt_cmdlineargparser_setint(&parser, "-i", 123, "some integer");
    assert(ier == 0);
    ier = mnt_cmdlineargparser_setdouble(&parser, "-d", 1.23, "some double");
    assert(ier == 0);
    ier = mnt_cmdlineargparser_setstring(&parser, "-s", "hello", "some string");
    assert(ier == 0);

    // parse the command line. Note: we need to pass a string array where
    // each element is a character array without the termination character

    const std::size_t string_len = 32; // max argument length
    char* args = new char[argc * string_len];
    for (size_t i = 0; i < (size_t) argc; ++i) {
        strncpy(&args[i * string_len], argv[i], string_len);
    }

    // parse 

    ier = mnt_cmdlineargparser_parse(&parser, argc, string_len, args);
    assert(ier == 0);

    // get the command argument values

    int ival;
    ier = mnt_cmdlineargparser_getint(&parser, "-i", &ival);
    assert(ier == 0);

    double dval;
    ier = mnt_cmdlineargparser_getdouble(&parser, "-d", &dval);
    assert(ier == 0);

    char sval[128];
    ier = mnt_cmdlineargparser_getstring(&parser, "-s", 128, sval);
    assert(ier == 0);

    // print results
    std::cout << "-i arg is " << ival << '\n';
    std::cout << "-d arg is " << dval << '\n';
    std::cout << "-s arg is " << sval << '\n';

    // clean up 

    ier = mnt_cmdlineargparser_del(&parser);
    assert(ier == 0);

    delete[] args;

    return 0;
}