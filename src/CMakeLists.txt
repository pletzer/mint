include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${NETCDF_INCLUDE_DIR}
)

set(LIB_FILES
  CmdLineArgParser.cpp
  MvVector.cpp
  MvMat_double.cpp
  mntLatLon.cpp
  mntPolysegmentIter.cpp
)

set(HEADER_FILES
  CmdLineArgParser.h
  mntLatLon.h
  mntGrid.h
  mntLineLineIntersector.h
  mntPolysegmentIter.h
)

add_library(mint ${LIB_FILES})
target_link_libraries(mint ${NETCDF_LIBRARIES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install library
install(TARGETS mint DESTINATION lib)
