###################################################
#
# CMakeLists.txt for MINT
#
###################################################

# Project information
project(mint)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})


# Required version
cmake_minimum_required(VERSION 3.1)

#
# Set default build type to Release
#
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE not specified -- setting to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)


set(DEBUG_CMAKE TRUE)
#
# Find VTK
#
if (VTK_INCLUDE_DIR AND VTK_LIBRARY_DIR AND VTK_VERSION)
    # allow the user to set the VTK libraries and include file location manually
    set(VTK_INCLUDE_DIR "" CACHE PATH "The directory containing the VTK include files")
    set(VTK_LIBRARY_DIR "" CACHE PATH "The directory containing the VTK library files")
    set(VTK_VERSION "" CACHE STRING "The VTK version, will search for libraries vtkIOCore-VTK_VERSION")
    set(VTK_LIBRARIES
        vtkIOCore-${VTK_VERSION}
        vtkIOLegacy-${VTK_VERSION}
        vtkCommonComputationalGeometry-${VTK_VERSION}
        vtkCommonExecutionModel-${VTK_VERSION}
        vtkCommonDataModel-${VTK_VERSION}
        vtkCommonTransforms-${VTK_VERSION}
        vtkCommonMisc-${VTK_VERSION}
        vtkCommonMath-${VTK_VERSION}
        vtkCommonSystem-${VTK_VERSION}
        vtkCommonCore-${VTK_VERSION}
        vtksys-${VTK_VERSION}
        )
    include_directories(${VTK_INCLUDE_DIR})
    link_directories(${VTK_LIBRARY_DIR})
    message(STATUS "User specified VTK location")
    message(STATUS "VTK include dir: ${VTK_INCLUDE_DIR}")
    message(STATUS "VTK library dir: ${VTK_LIBRARY_DIR}")
    message(STATUS "VTK version: ${VTK_VERSION}")
else ()
    # let MINT search for the location of VTK
    set(VTK_DIR "" CACHE PATH "The directory containing VTKConfig.cmake")
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    message(STATUS "VTK include dirs: ${VTK_INCLUDE_DIRS}")
    message(STATUS "VTK library dirs: ${VTK_LIBRARY_DIRS}")
endif()

# allow the user to set the variables manually
set(NETCDF_INCLUDE_DIR "" CACHE PATH "Directory containing the netCDF include files")
set(NETCDF_LIBRARIES "" CACHE STRING "List of netCDF libraries")

# Find the NetCDF library. Setting NETCDF_INCLUDE_DIR and NETCDF_LIBRARIES takes precedence.
# If not set then key off from nc-config
if (NETCDF_INCLUDE_DIR STREQUAL "" OR NETCDF_LIBRARIES STREQUAL "")
    execute_process(COMMAND "nc-config" "--includedir"
                    OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)
    execute_process(COMMAND "nc-config" "--libs"
                    OUTPUT_VARIABLE NETCDF_LIBRARIES)
endif()

if (NETCDF_INCLUDE_DIR STREQUAL "" OR NETCDF_LIBRARIES STREQUAL "")
    # we weree not able to find NetCDF
endif()

# remove \n at the end of the above variables
string(STRIP "${NETCDF_INCLUDE_DIR}" NETCDF_INCLUDE_DIR)
string(STRIP "${NETCDF_LIBRARIES}" NETCDF_LIBRARIES)


if ("${NETCDF_INCLUDE_DIR}" STREQUAL "" OR "${NETCDF_LIBRARIES}" STREQUAL "")
    message(STATUS "*******************************************************")
    message(STATUS "Unable to infer the location of NetCDF!")
    message(STATUS "Either:")
    message(STATUS "1) Specify -DNETCDF_INCLUDE_DIR:PATH=<dir>")
    message(STATUS "           -DNETCDF_LIBRARIES=\"-L<dir> -lnetcdf\"")
    message(STATUS "   (adding -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl and other libraries if needed)")
    message(STATUS "--or--")
    message(STATUS "2) Make sure to have nc-config in your PATH")
    message(STATUS "   so that NETCDF_INCLUDE_DIR and NETCDF_LIBRARIES can be inferred")
    message(FATAL_ERROR "*******************************************************")
endif()

include_directories(
  ${NETCDF_INCLUDE_DIR}
)

message(STATUS "NetCDF include dir: ${NETCDF_INCLUDE_DIR}")
message(STATUS "NetCDF libraries: ${NETCDF_LIBRARIES}")

find_package(Doxygen)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)



