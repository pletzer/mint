#ifdef _WIN32
#define LIBRARY_API extern "C"  __declspec(dllexport)
#else
#define LIBRARY_API extern "C"
#endif

