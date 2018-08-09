#include <mntCmdLineArgParser.h>
#include <cassert>
#include <cstring>
#include <iostream>

int main(int argc, char** argv) {

    CmdLineArgParser* self;
    int ier;

    // create 

    ier = mnt_cmdlineargparser_new(&self);
    assert(ier == 0);

    // set the command line options

    ier = mnt_cmdlineargparser_setint(&self, "-i", 123, "some integer");
    assert(ier == 0);
    ier = mnt_cmdlineargparser_setdouble(&self, "-d", 1.23, "some double");
    assert(ier == 0);
    ier = mnt_cmdlineargparser_setstring(&self, "-s", "hello", "some string");
    assert(ier == 0);

    // parse the command line. Note: we need too pass a string array where
    // each element is a character array without the termination character

    int* args_lengths = new int[argc];
    char** args = new char*[argc];
    for (size_t i = 0; i < (size_t) argc; ++i) {
        size_t n = strlen(argv[i]);
        args[i] = new char[n];
        strncpy(args[i], argv[i], n);
        args_lengths[i] = n;
    }

    // parse 

    ier = mnt_cmdlineargparser_parse(&self, argc, args_lengths, args);
    assert(ier == 0);

    // get the command argument values

    int ival;
    ier = mnt_cmdlineargparser_getint(&self, "-i", &ival);
    assert(ier == 0);

    double dval;
    ier = mnt_cmdlineargparser_getdouble(&self, "-d", &dval);
    assert(ier == 0);

    char sval[128];
    int n;
    ier = mnt_cmdlineargparser_getstring(&self, "-s", sval, &n);
    assert(ier == 0);

    // print results
    std::cout << "-i arg is " << ival << '\n';
    std::cout << "-d arg is " << dval << '\n';
    std::cout << "-s arg is " << sval << '\n';

    // clean up 

    ier = mnt_cmdlineargparser_del(&self);
    assert(ier == 0);

    for (size_t i = 0; i < (size_t) argc; ++i) {
        delete[] args[i];
    }
    delete[] args;
    delete[] args_lengths;

    return 0;
}