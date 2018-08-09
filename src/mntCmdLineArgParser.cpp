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
int mnt_cmdlineargparser_parse(CmdLineArgParser** self, int nargs, const int args_lengths[], char** args) {
	char** argv = new char*[nargs];
	for (size_t i = 0; i < (size_t) nargs; ++i) {
		argv[i] = new char[args_lengths[i]];
	}
	(*self)->parse(nargs, argv);
	for (size_t i = 0; i < (size_t) nargs; ++i) {
		delete[] argv[i];
	}
	delete[] argv;
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
int mnt_cmdlineargparser_getstring(CmdLineArgParser** self, const char* name, char* val, int* n){
	size_t nv = strlen(val);
	std::cerr << "** nv = " << nv << '\n';
	std::string sval = (*self)->get<std::string>(name);
	*n = (int) sval.size();
	strncpy(val, sval.c_str(), nv);
	std::cerr << "val=" << val << "< sval=" << sval << '\n';
	if (*n > nv) {
		return 1;
	}
	return 0;
}
