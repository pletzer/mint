include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${NETCDF_INCLUDE_DIR}
)

set(LIB_FILES
  mntUgrid2D.cpp
  mntQuadEdgeIter.cpp
  mntPolylineParser.cpp
  mntLineGridIntersector.cpp
  CmdLineArgParser.cpp
  GrExprParser.cpp
  GrExprAdaptor.cpp
  MvVector.cpp
  MvMat_double.cpp
  mntLatLon.cpp
  mntGrid.cpp
  mntPolysegmentIter.cpp
  mntLineTriangleIntersector.cpp
  mntRegridEdges.cpp
  mntRegridEdges2D.cpp
  mntCellLocator.cpp
  mntCmdLineArgParser.cpp
  mnt_celllocator_capi_mod.f90
  mnt_regridedges_capi_mod.f90
  mnt_cmdlineargparser_capi_mod.f90
)

set(HEADER_FILES
  mntUgrid2D.h
  mntQuadEdgeIter.h
  mntPolylineParser.h
  mntLineGridIntersector.h
  CmdLineArgParser.h
  GrExprParser.h
  GrExprAdaptor.h
  mntLatLon.h
  mntGrid.h
  mntLineLineIntersector.h
  mntPolysegmentIter.h
  mntRegridEdges.h
  mntRegridEdges2D.h
  mntCellLocator.h
  mntCmdLineArgParser.h
  MvFunctors.h MvMatrix.h MvVector.h
)

add_library(mint ${LIB_FILES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ".")

# Install library
install(TARGETS mint DESTINATION lib)
