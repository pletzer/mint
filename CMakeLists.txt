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
set(CMAKE_CXX_STANDARD 11)

set(DEBUG_CMAKE TRUE)

find_package(VTK REQUIRED)
#find_package(NetCDF REQUIRED)
# the above does not work well,
# user must explicitly specify 
# NETCDF_DIR and we will assume 
# $NETCDF_DIR/include and NETCDF_DIR/lib
# to be below
find_package(HDF5)
find_package(CURL)

if (ZLIB_FOUND)
	message(STATUS "ZLIB_LIBRARIES are ${ZLIB_LIBRARIES}")
endif()
if (NetCDF_FOUND)
	message(STATUS "NetCDF libraries are >${netCDF_LIBRARIES}< >${NETCDF_LIBRARY}< >${NETCDF_LIBRARIES}<")
endif()

include(${VTK_USE_FILE})

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)
