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

# Add version to each library name
set(VTK_VER "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
#foreach (LIB in ${VTK_LIBRARIES})
#	string(REGEX REPLACE "(vtk[^ ]+)" "\\1-${VTK_VER}" LIB2 "${LIB}")
#	message(STATUS "${LIB} ${LIB2}")
#endforeach()

# Replace ; with space
string(REPLACE ";" " " VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES}")

# remove full path libraries, containing dylib, so, etc.
string(REGEX REPLACE "([^ ]+\.dylib)" "" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")
#string(REGEX REPLACE "([^ ]+\.so)" "" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")
#string(REGEX REPLACE "([^ ]+\.a)" "" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")

# Remove non vtk libraries
string(REGEX REPLACE "verdict" "" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")

# remove Python libraries
string(REGEX REPLACE "([^ ]+Python[^ ]+)" "" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")
# Add version to library names
string(REGEX REPLACE "(vtk[^ ]+)" "\\1-${VTK_VER}" VTK_LIBRARIES_VERSION_SPACE "${VTK_LIBRARIES_VERSION_SPACE}")
message(STATUS "${VTK_LIBRARIES_VERSION_SPACE}")

# Create file containing list of libraries, save in mint.mk
configure_file(mint.mk.in mint.mk)
install(FILES "${CMAKE_BINARY_DIR}/mint.mk" DESTINATION share)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)
