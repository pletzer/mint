###################################################
#
# CMakeLists.txt for MINT
#
###################################################

# Project information
project(mint)

# Required version
cmake_minimum_required(VERSION 3.5)

# get the MINT version
file(READ "${CMAKE_SOURCE_DIR}/version.txt" VERSION)

#
# Set default build type to Release
#
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE not specified -- setting to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)
enable_testing()

option(BUILD_FORTRAN "Whether or not to build Fortran code" OFF)
if (BUILD_FORTRAN)
    enable_language(Fortran)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

set(DEBUG_CMAKE TRUE)
#
# Find VTK
#
# allow the user to set the variables manually
set(VTK_INCLUDE_DIR "" CACHE PATH "Directories containing the VTK header")
set(VTK_LIBRARY_DIR "" CACHE PATH "Directories containing the VTK libraries")
set(VTK_VERSION "" CACHE STRING "The VTK version, will search for libraries vtkXXX-VTK_VERSION")
if ((NOT "${VTK_INCLUDE_DIR}" STREQUAL "") AND (NOT "${VTK_LIBRARY_DIR}" STREQUAL ""))

    # infer the VTK version from the library names, expecting *vtkIOCore-VTK_VERSION.*
    # fetch a library file
    file(GLOB_RECURSE f "${VTK_LIBRARY_DIR}/*vtkIOCore-[0-9].[0-9].[a-zA-Z]*")
    # get rid of the directory
    get_filename_component(bf "${f}" NAME)
    # search for the version string, X.Y in libvtkIOCore-X.Y.dylib
    string(REGEX MATCH "[0-9]\.[0-9]+" VTK_VERSION "${bf}")
    message(STATUS "VTK_VERSION = ${VTK_VERSION}")

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
    set(VTK_DIR "" CACHE PATH "The directory containing VTKConfig.cmake or vtk-config.cmake")
    if (VTK_VERSION VERSION_LESS "8.90.0")
        find_package(VTK REQUIRED)
    else()
        message(STATUS "Using new VTK module system")
        find_package(VTK REQUIRED COMPONENTS
            vtksys
            CommonCore
            CommonSystem
            CommonMath
            CommonMisc
            CommonTransforms
            CommonDataModel
            CommonExecutionModel
            CommonComputationalGeometry
            IOLegacy
            IOCore
        )
    endif()
endif()


#
# Find NetCDF
#

# allow the user to set the variables manually
set(NETCDF_INCLUDE_DIRS "" CACHE PATH "Directories containing the netCDF C include files and Fortran module files")
set(NETCDF_LIBRARIES "" CACHE STRING "List of netCDF libraries")

# Setting NETCDF_INCLUDE_DIRS and NETCDF_LIBRARIES takes precedence.
# If not set then key off from nc-config and optionally nf-config
if (NETCDF_INCLUDE_DIRS STREQUAL "" OR NETCDF_LIBRARIES STREQUAL "")
    execute_process(COMMAND "nc-config" "--includedir"
        OUTPUT_VARIABLE C_NETCDF_INCLUDE_DIR)
    execute_process(COMMAND "nc-config" "--libs"
        OUTPUT_VARIABLE C_NETCDF_LIBRARIES)
    # remove \n at the end of the above variables
    string(STRIP "${C_NETCDF_INCLUDE_DIR}" C_NETCDF_INCLUDE_DIR)
    string(STRIP "${C_NETCDF_LIBRARIES}" C_NETCDF_LIBRARIES)

    if (BUILD_FORTRAN)
        execute_process(COMMAND "nf-config" "--includedir"
            OUTPUT_VARIABLE FORTRAN_NETCDF_INCLUDE_DIR)
        execute_process(COMMAND "nf-config" "--flibs"
            OUTPUT_VARIABLE FORTRAN_NETCDF_LIBRARIES)
        if (FORTRAN_NETCDF_INCLUDE_DIR STREQUAL "" OR FORTRAN_NETCDF_LIBRARIES STREQUAL "")
            message(STATUS "Command nf-config was not found")
            set(BUILD_FORTRAN OFF)
        else()
            string(STRIP "${FORTRAN_NETCDF_INCLUDE_DIR}" FORTRAN_NETCDF_INCLUDE_DIR)
            string(STRIP "${FORTRAN_NETCDF_LIBRARIES}" FORTRAN_NETCDF_LIBRARIES)
        endif()
    endif()

    set(NETCDF_INCLUDE_DIRS ${C_NETCDF_INCLUDE_DIR} ${FORTRAN_NETCDF_INCLUDE_DIR})
    set(NETCDF_LIBRARIES ${C_NETCDF_LIBRARIES} ${FORTRAN_NETCDF_LIBRARIES})
    if (NOT BUILD_FORTRAN)
        message(STATUS "Will not build Fortran interface")
    endif()
endif()

if ("${NETCDF_INCLUDE_DIRS}" STREQUAL "" OR "${NETCDF_LIBRARIES}" STREQUAL "")
    message(STATUS "*******************************************************")
    message(STATUS "Unable to infer the location of NetCDF!")
    message(STATUS "Either:")
    message(STATUS "1) Specify -DNETCDF_INCLUDE_DIRS:PATH=<dir>")
    message(STATUS "           -DNETCDF_LIBRARIES=\"-L<dir> [-lnetcdff] -lnetcdf\"")
    message(STATUS "   (adding -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl and other libraries if needed)")
    message(STATUS "--or--")
    message(STATUS "2) Make sure to have nc-config [and optionally nf-config] in your PATH")
    message(STATUS "   so that NETCDF_INCLUDE_DIRS and NETCDF_LIBRARIES can be inferred")
    message(FATAL_ERROR "*******************************************************")
endif()

include_directories(
    ${NETCDF_INCLUDE_DIRS}
)

message(STATUS "NetCDF include dirs: ${NETCDF_INCLUDE_DIRS}")
message(STATUS "NetCDF libraries: ${NETCDF_LIBRARIES}")

#
# Doxygen (optional)
#
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxygen/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS}")
if (BUILD_FORTRAN)
    message(STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
endif()
