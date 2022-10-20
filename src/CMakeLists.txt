include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${NETCDF_INCLUDE_DIR}
)

set(LIB_FILES
  mntLogger.cpp
  mntFileMeshNameExtractor.cpp
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
  mntPolylineIntegral.cpp
  mntWeights.cpp
  mntVectorInterp.cpp
  mntExtensiveFieldConverter.cpp
  mntExtensiveFieldAdaptor.cpp
)

if (BUILD_FORTRAN)
  set(LIB_FILES_FORTRAN
    mnt_logger_capi_mod.f90
    mnt_ncfieldread_capi_mod.f90
    mnt_ncfieldwrite_capi_mod.f90
    mnt_celllocator_capi_mod.f90
    mnt_regridedges_capi_mod.f90
    mnt_cmdlineargparser_capi_mod.f90
    mnt_multiarrayiter_capi_mod.f90
    mnt_ncdimensions_capi_mod.f90
    mnt_ncattributes_capi_mod.f90
 )
endif()

set(HEADER_FILES
  mntLIBRARY_API.h
  mntGlobal.h
  mntLogger.h
  mntFileMeshNameExtractor.h
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
  mntPolylineIntegral.h
  mntMatMxN.h
  MvFunctors.h MvMatrix.h MvVector.h
  mntVecN.h
  mntWeights.h
  mntVectorInterp.h
  mntExtensiveFieldConverter.h
  mntExtensiveFieldAdaptor.h
)

if (VTK_VERSION VERSION_LESS "8.90.0")
    include(${VTK_USE_FILE})
    add_library(mint ${LIB_FILES} ${LIB_FILES_FORTRAN})
else()
    message(STATUS "Using new VTK module system")
    add_library(mint ${LIB_FILES} ${LIB_FILES_FORTRAN})
    # may want to replace ${VTK_LIBRARIES} with VTK::vtksys VTK::CommonCore etc
    target_link_libraries(mint ${visibility} ${VTK_LIBRARIES})
    #vtk_module_autoinit(TARGETS mint MODULES ${VTK_LIBRARIES})
endif()

# Install headers
install(FILES ${HEADER_FILES} DESTINATION include)

# Install fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ".")

# Install library
install(TARGETS mint DESTINATION lib)
