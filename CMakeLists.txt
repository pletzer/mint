###################################################
#
# CMakeLists.txt for MINT
#
###################################################

# Project information
project(mint)
set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 3.1)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

set(DEBUG_CMAKE TRUE)

#
# Find VTK
#
# allow the user to set the variables manually
set(VTK_DIR "" CACHE PATH "The directory containing VTKConfig.cmake")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
message(STATUS "VTK include dirs: ${VTK_INCLUDE_DIRS}")
message(STATUS "VTK library dirs: ${VTK_LIBRARY_DIRS}")


#
# Find the netcdf library and its dependencies by keying off from the nc-config
# command. 
#
set(NETCDF_LIBRARY_DIR "" CACHE PATH "The top directory where NetCDF is installed")

# allow the user to set the variables manually
set(NETCDF_INCLUDE_DIR "" CACHE PATH "Directory containing the netCDF include files")
set(NETCDF_LIBRARIES "" CACHE PATH "List of netCDF libraries")

if (NETCDF_INCLUDE_DIR STREQUAL "")
    execute_process(COMMAND "nf-config" "--includedir"
                    OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)
endif()

if (NETCDF_LIBRARIES STREQUAL "")
    string(STRIP "${NETCDF_Fortran_LIBRARIES}" NETCDF_Fortran_LIBRARIES)
    execute_process(COMMAND "nc-config" "--libs"
                    OUTPUT_VARIABLE NETCDF_C_LIBRARIES)
    string(STRIP "${NETCDF_C_LIBRARIES}" NETCDF_C_LIBRARIES)
    set(NETCDF_LIBRARIES "${NETCDF_C_LIBRARIES}")
endif()

if (NOT NETCDF_LIBRARY_DIR STREQUAL "")
    # sometimes the netcdf library path returned by nc-config is wrong
    # so users should be able to set expicitly
    set(NETCDF_LIBRARIES "-L ${NETCDF_LIBRARY_DIR} ${NETCDF_LIBRARIES}")
endif()

# remove \n at the end of the above variables
string(STRIP "${NETCDF_INCLUDE_DIR}" NETCDF_INCLUDE_DIR)
string(STRIP "${NETCDF_LIBRARIES}" NETCDF_LIBRARIES)


if ("${NETCDF_INCLUDE_DIR}" STREQUAL "" OR "${NETCDF_LIBRARIES}" STREQUAL "")
    message(STATUS "*******************************************************")
    message(STATUS "Unable to infer the location of NetCDF!")
    message(STATUS "Either:")
    message(STATUS "1) Specify -DNETCDF_INCLUDE_DIR:PATH=<dir>")
    message(STATUS "           -DNETCDF_LIBRARIES:PATH=\"-L<dir> -lnetcdff -lnetcdf\"")
    message(STATUS "   (adding -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl and other libraries if need be)")
    message(STATUS "--or--")
    message(STATUS "2) Make sure to have nf-config and nc-config in your PATH")
    message(STATUS "   so that NETCDF_INCLUDE_DIR and NETCDF_LIBRARIES can be deduced")
    message(FATAL_ERROR "*******************************************************")
endif()

include_directories(
  ${NETCDF_INCLUDE_DIR}
)


enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)


