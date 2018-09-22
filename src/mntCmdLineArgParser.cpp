#include <mntCmdLineArgParser.h>
#include <cstring>

extern "C"
int mnt_cmdlineargparser_new(CmdLineArgParser** self) {
	*self = new CmdLineArgParser();
	return 0;
}

extern "C"
int mnt_cmdlineargparser_del(CmdLineArgParser** self) {
	delete *self;
	return 0;
}

extern "C"
int mnt_cmdlineargparser_setdouble(CmdLineArgParser** self, const char* name, double def_value, const char* help) {
	(*self)->set(name, def_value, help);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_setint(CmdLineArgParser** self, const char* name, int def_value, const char* help) {
	(*self)->set(name, def_value, help);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_setstring(CmdLineArgParser** self, const char* name, const char* def_value, const char* help) {
	(*self)->set(name, std::string(def_value), help);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_setbool(CmdLineArgParser** self, const char* name, int def_value, const char* help) {
	bool v = false;
	if (def_value != 0) v = true;
	(*self)->set(name, v, help);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_parse(CmdLineArgParser** self, int nargs, int n, char* args) {

	char** argv = new char*[nargs];
	for (size_t i = 0; i < (size_t) nargs; ++i) {
		argv[i] = new char[n];
		// assumes args is contiguous in memorygs 
		strncpy(argv[i], &args[n*i], n);
	}

	bool success = (*self)->parse(nargs, argv);

	// clean up
	for (size_t i = 0; i < (size_t) nargs; ++i) {
		delete[] argv[i];
	}
	delete[] argv;

	if (!success) return 1;
	return 0;
}

extern "C"
int mnt_cmdlineargparser_help(CmdLineArgParser** self) {
	(*self)->help();
	return 0;
}


extern "C"
int mnt_cmdlineargparser_getdouble(CmdLineArgParser** self, const char* name, double* val) {
	*val = (*self)->get<double>(name);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_getint(CmdLineArgParser** self, const char* name, int* val){
	*val = (*self)->get<int>(name);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_getstring(CmdLineArgParser** self, const char* name, int n, char* val){
	std::string sval = (*self)->get<std::string>(name);
	if (sval.size() == 0) {
		// not valid name
		return 1;
	}
    // initialize the receiving string to the termination character
    for (int i = 0; i < n; ++i) {
        val[i] = '\0';
    }
    // copy
	int n2 = std::min((int) sval.size(), n);
	strncpy(val, sval.c_str(), n2);
	return 0;
}

extern "C"
int mnt_cmdlineargparser_getbool(CmdLineArgParser** self, const char* name, int* val) {
	bool v = (*self)->get<bool>(name);
	if (v) {
		*val = 1;
	}
	else {
		*val = 0;
	}
	return 0;
}

