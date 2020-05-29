include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${NETCDF_INCLUDE_DIR}
)

set(LIB_FILES
  vmtCellLocator.cpp
  mntUgrid2D.cpp
  mntQuadEdgeIter.cpp
  mntPolylineParser.cpp
  CmdLineArgParser.cpp
  GrExprParser.cpp
  GrExprAdaptor.cpp
  MvVector.cpp
  mntVecN.cpp
  MvMat_double.cpp
  mntMatMxN.cpp
  mntLatLon.cpp
  mntGrid.cpp
  mntPolysegmentIter.cpp
  mntLineTriangleIntersector.cpp
  mntRegridEdges.cpp
  mntCellLocator.cpp
  mntCmdLineArgParser.cpp
  mntNcAttributes.cpp
  mntNcFieldRead.cpp
  mntNcFieldWrite.cpp
  mntNcDimensions.cpp
  mntMultiArrayIter.cpp
  mntRegridAxis.cpp
  mnt_ncfieldread_capi_mod.f90
  mnt_ncfieldwrite_capi_mod.f90
  mnt_celllocator_capi_mod.f90
  mnt_regridedges_capi_mod.f90
  mnt_cmdlineargparser_capi_mod.f90
)

set(HEADER_FILES
  vmtCellLocator.h
  mntUgrid2D.h
  mntQuadEdgeIter.h
  mntPolylineParser.h
  CmdLineArgParser.h
  GrExprParser.h
  GrExprAdaptor.h
  mntLatLon.h
  mntGrid.h
  mntLineLineIntersector.h
  mntPolysegmentIter.h
  mntRegridEdges.h
  mntCellLocator.h
  mntCmdLineArgParser.h
  mntNcAttributes.h
  mntNcFieldRead.h
  mntNcFieldWrite.h
  mntNcDimensions.h
  mntMultiArrayIter.h
  mntRegridAxis.h
  mntMatMxN.h
  MvFunctors.h MvMatrix.h MvVector.h
  mntVecN.h
)

add_library(mint ${LIB_FILES})

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ".")

# Install library
install(TARGETS mint DESTINATION lib)
