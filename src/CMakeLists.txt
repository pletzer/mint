include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${NETCDF_INCLUDE_DIR}
)

set(LIB_FILES
  vmtCellLocator.cpp
  mntUgrid2D.cpp
  mntQuadEdgeIter.cpp
  mntPolylineParser.cpp
  mntLineGridIntersector.cpp
  CmdLineArgParser.cpp
  GrExprParser.cpp
  GrExprAdaptor.cpp
  MvVector.cpp
  mntVec2.cpp mntVec3.cpp mntVec4.cpp mntVec9.cpp
  MvMat_double.cpp
  mntMat3x3.cpp mntMat2x2.cpp
  mntLatLon.cpp
  mntGrid.cpp
  mntPolysegmentIter.cpp
  mntLineTriangleIntersector.cpp
  mntRegridEdges.cpp
  mntRegridEdges2D.cpp
  mntCmdLineArgParser.cpp
  mnt_celllocator_capi_mod.f90
  mnt_regridedges_capi_mod.f90
  mnt_cmdlineargparser_capi_mod.f90
)

set(HEADER_FILES
  vmtCellLocator.h
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
  mntCmdLineArgParser.h
  mntMat3x3.h mntMat2x2.h
  MvFunctors.h MvMatrix.h MvVector.h
  mntVec2.h mntVec3.h mntVec4.h mntVec9.h
)

add_library(mint ${LIB_FILES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ".")

# Install library
install(TARGETS mint DESTINATION lib)
