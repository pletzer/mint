include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

configure_file(testSimpleRegridEdges.cpp testSimpleRegridEdges.cxx)
configure_file(testUgridLocator.cpp testUgridLocator.cxx)
configure_file(testNetcdf.cpp testNetcdf.cxx)
configure_file(testNcAttributes.cpp testNcAttributes.cxx)
configure_file(testNcDimensions.cpp testNcDimensions.cxx)
configure_file(testUgrid2D.cpp testUgrid2D.cxx)
configure_file(testRegridEdgesFromUgrid.cpp testRegridEdgesFromUgrid.cxx)
configure_file(testGrid.cpp testGrid.cxx)
configure_file(testRegridEdgesFromUgrid.cpp testRegridEdgesFromUgrid.cxx)
configure_file(test_ncfieldread.f90.in test_ncfieldread.f90)
configure_file(test_ncdimensions.f90.in test_ncdimensions.f90)


add_executable(test_ncdimensions ${CMAKE_BINARY_DIR}/tests/test_ncdimensions.f90)
set_target_properties(test_ncdimensions PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_ncdimensions
                      mint
                      ${VTK_LIBRARIES}
                      netcdff ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(test_multiarrayiterF test_multiarrayiterF.f90)
set_target_properties(test_multiarrayiterF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_multiarrayiterF
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)


add_executable(test_ncfieldwrite test_ncfieldwrite.f90)
set_target_properties(test_ncfieldwrite PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_ncfieldwrite
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(test_ncfieldread ${CMAKE_BINARY_DIR}/tests/test_ncfieldread.f90)
set_target_properties(test_ncfieldread PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_ncfieldread
                      mint
                      ${VTK_LIBRARIES}
                      netcdff ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(test_cell_locator_capi test_cell_locator_capi.f90)
set_target_properties(test_cell_locator_capi PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(test_cell_locator_capi
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testVecMat testVecMat.cxx)
target_link_libraries(testVecMat
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testUgridLocator testUgridLocator.cxx)
target_link_libraries(testUgridLocator
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testNetcdf testNetcdf.cxx)
target_link_libraries(testNetcdf
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testNcAttributes testNcAttributes.cxx)
target_link_libraries(testNcAttributes
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testNcDimensions testNcDimensions.cxx)
target_link_libraries(testNcDimensions
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testMultiArrayIter testMultiArrayIter.cxx)
target_link_libraries(testMultiArrayIter
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testUgrid2D testUgrid2D.cxx)
target_link_libraries(testUgrid2D
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testSimpleRegridEdges testSimpleRegridEdges.cxx)
target_link_libraries(testSimpleRegridEdges
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testRegridEdgesFromUgrid testRegridEdgesFromUgrid.cxx)
target_link_libraries(testRegridEdgesFromUgrid
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testPolylineParser testPolylineParser.cxx)
target_link_libraries(testPolylineParser
                    mint
                    ${VTK_LIBRARIES}
)

add_executable(testLineLineIntersector testLineLineIntersector.cxx)
target_link_libraries(testLineLineIntersector
                    mint
                    ${VTK_LIBRARIES}
)

add_executable(testLineTriangleIntersector testLineTriangleIntersector.cxx)
target_link_libraries(testLineTriangleIntersector
                    mint
                    ${VTK_LIBRARIES}
)

add_executable(testPolysegmentIter testPolysegmentIter.cxx)
target_link_libraries(testPolysegmentIter
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testGrid testGrid.cxx)
target_link_libraries(testGrid
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testCellLocator testCellLocator.cxx)
target_link_libraries(testCellLocator
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testC2fString testC2fString.f90)
set_target_properties(testC2fString PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testC2fString
                      mint
                      ${VTK_LIBRARIES}
)


add_executable(testExprParser testExprParser.cxx)
target_link_libraries(testExprParser
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testSize testSize.cxx)
target_link_libraries(testSize
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testCmdLineArgParserC testCmdLineArgParserC.cxx)
target_link_libraries(testCmdLineArgParserC
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testRegridEdgesF testRegridEdgesF.f90)
set_target_properties(testRegridEdgesF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testRegridEdgesF
                    mint
                    ${VTK_LIBRARIES}
                    ${NETCDF_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${CURL_LIBRARIES}
)

add_executable(testCmdLineArgParserF testCmdLineArgParserF.f90)
set_target_properties(testCmdLineArgParserF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testCmdLineArgParserF
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testCmdLineArgParserFB testCmdLineArgParserFB.f90)
set_target_properties(testCmdLineArgParserFB PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testCmdLineArgParserFB
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(testRegridAxis testRegridAxis.cxx)
target_link_libraries(testRegridAxis
                      mint
                      ${VTK_LIBRARIES}
)

add_test(NAME test_multiarrayiterF
         COMMAND test_multiarrayiterF)
set_tests_properties(test_multiarrayiterF PROPERTIES
                     FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME test_ncdimensionsF
         COMMAND test_ncdimensions)
set_tests_properties(test_ncdimensionsF PROPERTIES
                     FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME test_ncfieldwrite
         COMMAND test_ncfieldwrite)
set_tests_properties(test_ncfieldwrite PROPERTIES
                     FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME test_ncfieldread
         COMMAND test_ncfieldread)
set_tests_properties(test_ncfieldread PROPERTIES
                     FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME test_cell_locator_capi
         COMMAND test_cell_locator_capi)

add_test(NAME regridEdgesF
         COMMAND testRegridEdgesF 
         "-s" "${CMAKE_SOURCE_DIR}/data/cs_16.nc:physics" 
         "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics" 
         "-o" "weights.nc")
set_tests_properties(regridEdgesF PROPERTIES
                     FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME vecMat COMMAND testVecMat)
add_test(NAME ugridLocator COMMAND testUgridLocator)
add_test(NAME netcdf COMMAND testNetcdf)
add_test(NAME attributes COMMAND testNcAttributes)
add_test(NAME dimensions COMMAND testNcDimensions)
add_test(NAME multiarrayiter COMMAND testMultiArrayIter)
add_test(NAME ugrid2D_1 COMMAND testUgrid2D 
                        "-s" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
                        "-vtk" "cs_4.vtk")
add_test(NAME ugrid2D_2 COMMAND testUgrid2D 
                        "-s" "${CMAKE_SOURCE_DIR}/data/output_UGRID.nc:Mesh2d"
                        "-vtk" "output_UGRID.vtk")
add_test(NAME simpleRegridEdges COMMAND testSimpleRegridEdges)
add_test(NAME regridEdgesFromUgrid COMMAND testRegridEdgesFromUgrid)
add_test(NAME polylineParser COMMAND testPolylineParser)
add_test(NAME lineLineIntersector COMMAND testLineLineIntersector)
add_test(NAME lineTriangleIntersector COMMAND testLineTriangleIntersector)
add_test(NAME polysegmentIter COMMAND testPolysegmentIter)
add_test(NAME grid COMMAND testGrid)
add_test(NAME cellLocator COMMAND testCellLocator)
add_test(NAME exprParser COMMAND testExprParser)
add_test(NAME size COMMAND testSize)
add_test(NAME cmdLineArgParserC COMMAND testCmdLineArgParserC)
add_test(NAME cmdLineArgParserC2 
        COMMAND testCmdLineArgParserC "-i" "-234" "-d" -2.34 "-s" "hi_there")
add_test(NAME cmdLineArgParserF COMMAND testCmdLineArgParserF)
add_test(NAME cmdLineArgParserF2 
         COMMAND testCmdLineArgParserF "-i" "input.vtk" "-o" "output.vtk" "-d" "-2.34" "-n" "100" "-v")
add_test(NAME cmdLineArgParserFB
         COMMAND testCmdLineArgParserFB "-i" "in.vtk")
add_test(NAME c2fString
         COMMAND testC2fString)

add_test(NAME regridedges_edge_field_lfric_diag
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/lfric_diag_ex.nc:Mesh2d_half_levels" "-S" "1"
            "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics" "-D" "1"
            "-v" "edge_field@${CMAKE_SOURCE_DIR}/data/edge_field_lfric_diag_ex.nc"
            "-O" "regridedges_edge_field_lfric_diag_ex.nc"
            "-o" "regridedges_edge_field_lfric_diag_ex.vtk")
set_tests_properties(regridedges_edge_field_lfric_diag PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-0[8-9]")

add_test(NAME regridedges_xyzt
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/lonlatzt_8x4x3x2.nc:mesh2d" "-S" "0"
            "-d" "${CMAKE_SOURCE_DIR}/data/c24_u_integrated.nc:physics" "-D" "1"
            "-v" "u@${CMAKE_SOURCE_DIR}/data/lonlatzt_8x4x3x2.nc"
            "-O" "regridedges_xyzt.nc"
            "-o" "regridedges_xyzt.vtk")
set_tests_properties(regridedges_xyzt PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regridedges_xyzt2
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/lonlatzt_100x50x3x2.nc:mesh2d" "-S" "0"
            "-d" "${CMAKE_SOURCE_DIR}/data/c24_u_integrated.nc:physics" "-D" "1"
            "-v" "u@${CMAKE_SOURCE_DIR}/data/lonlatzt_100x50x3x2.nc"
            "-O" "regridedges_xyzt2.nc"
            "-o" "regridedges_xyzt2.vtk")
set_tests_properties(regridedges_xyzt2 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regridedges_axes
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/lfric_diag_ex.nc:Mesh2d_half_levels" "-S" "1"
            "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics" "-D" "1"
            "-v" "u1@${CMAKE_SOURCE_DIR}/data/lfric_diag_ex.nc"
            "-O" "regridedges_axes.nc"
            "-o" "regridedges_axes.vtk")
# we do not expect conseervation
#set_tests_properties(regridedges_axes PROPERTIES
#                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regridedges_0
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc:latlon" "-S" "0"
            "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics" "-D" "1"
            "-v" "edge_integrated_velocity@${CMAKE_SOURCE_DIR}/data/latlon4x2.nc"
            "-O" "regridedges_0_data.nc"
            "-o" "regridedges_0.vtk")
set_tests_properties(regridedges_0 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regridedges_1
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc:latlon" "-S" "0"
            "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics" "-D" "1"
            "-w" "regridedges_1_weights.nc")
set_tests_properties(regridedges_1 PROPERTIES
                     PASS_REGULAR_EXPRESSION "info: saving weights in file")

add_test(NAME regridedges_2
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
            "-s" "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc:latlon"
            "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
            "-W" "${CMAKE_SOURCE_DIR}/data/regridedges_1_weights.nc"
            "-v" "edge_integrated_velocity@${CMAKE_SOURCE_DIR}/data/latlon4x2.nc"
            "-O" "regridedges_2_data.nc"
            "-o" "regridedges_2.vtk")
set_tests_properties(regridedges_2 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesLATLON4x2ToCS4
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges"
                 "-S" "0"
                 "-D" "1"
                 "-N" "1"
                 "-s" "${CMAKE_SOURCE_DIR}/data/latlon4x2.nc:latlon"
                 "-v" "edge_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
                 "-o" "regrid_edges_output_LATLON4x2ToCS4.vtk")
set_tests_properties(regrid_edgesLATLON4x2ToCS4 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesLATLON100x50ToCS4_1cellPerBucket
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc:latlon"
                 "-S" "0" 
                 "-v" "edge_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
                 "-o" "regrid_edges_output_LATLON100x50ToCS4.vtk"
                 "-N" "1")
set_tests_properties(regrid_edgesLATLON100x50ToCS4_1cellPerBucket PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesLATLON100x50ToCS4
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc:latlon"
                 "-S" "0" 
                 "-v" "edge_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
                 "-o" "regrid_edges_output_LATLON100x50ToCS4.vtk")
set_tests_properties(regrid_edgesLATLON100x50ToCS4 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesLATLON100x50ToCS16
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-N" "20"
                 "-s" "${CMAKE_SOURCE_DIR}/data/latlon100x50.nc:latlon"
                 "-S" "0" 
                 "-v" "edge_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_16.nc:physics"
                 "-o" "regrid_edges_output_LATLON100x50ToCS16.vtk")
set_tests_properties(regrid_edgesLATLON100x50ToCS4 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesUgrid16To4
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/vel_cs_16.nc:physics" 
                 "-v" "line_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_4.nc:physics"
                 "-o" "regrid_edges_output_Ugrid16To4.vtk")
set_tests_properties(regrid_edgesUgrid16To4 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_edgesUgrid16To16
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/vel_cs_16.nc:physics" 
                 "-v" "line_integrated_velocity" 
                 "-d" "${CMAKE_SOURCE_DIR}/data/cs_16.nc:physics"
                 "-o" "regrid_edges_output_Ugrid16To4.vtk")
set_tests_properties(regrid_edgesUgrid16To4 PROPERTIES
                     PASS_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/[^/]*/[^e]+e-1[0-9]")

add_test(NAME regrid_nan
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/grid.nc:physics" "-S" "1"
                 "-d" "${CMAKE_SOURCE_DIR}/data/upstream.nc:physics" "-D" "1"
                 "-v" "edge_integrated_velocity@${CMAKE_SOURCE_DIR}/data/input_data_initial.nc"
                 "-o" "regrid_nan.vtk"
                 "-O" "regrid_nan_data.nc"
                 )
set_tests_properties(regrid_nan PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/nan/")

add_test(NAME regrid_nan_load_weights
         COMMAND "${CMAKE_BINARY_DIR}/tools/regrid_edges" 
                 "-s" "${CMAKE_SOURCE_DIR}/data/grid.nc:physics" "-S" "1"
                 "-d" "${CMAKE_SOURCE_DIR}/data/upstream.nc:physics" "-D" "1"
                 "-W" "${CMAKE_SOURCE_DIR}/data/weights.nc"
                 "-v" "edge_integrated_velocity@${CMAKE_SOURCE_DIR}/data/input_data_initial.nc"
                 "-o" "regrid_nan_load_weights.vtk"
                 "-O" "regrid_nan_load_weights_data.nc"
                 )
set_tests_properties(regrid_nan_load_weights PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Min/avg/max cell loop integrals: [^/]*/nan/")

add_test(NAME line_proj1
         COMMAND "${CMAKE_BINARY_DIR}/tools/line_proj" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/latlon4x2.vtk"
                 "-v" "edge_integrated_velocity"
                 "-N" "1"
                 "-p" "(0., 22.5), (22.5, 20.9410204722438422)"
                 "-verbose")
set_tests_properties(line_proj1 PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total integrated length for segment")

add_test(NAME line_proj_inside
         COMMAND "${CMAKE_BINARY_DIR}/tools/line_proj" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/latlon4x2.vtk"
                 "-v" "edge_integrated_velocity"
                 "-N" "1"
                 "-p" "(359., 22.5), (337.5, 20.9410204722438422)"
                 "-P" "360."
                 "-verbose")
set_tests_properties(line_proj_inside PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total integrated length for segment")

add_test(NAME line_proj_partially_outside
         COMMAND "${CMAKE_BINARY_DIR}/tools/line_proj" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/latlon4x2.vtk"
                 "-v" "edge_integrated_velocity"
                 "-N" "1"
                 "-p" "(-1., 22.5), (22.5, 20.9410204722438422)"
                 "-P" "360."
                 "-verbose")
set_tests_properties(line_proj_partially_outside PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total integrated length for segment")

add_test(NAME line_proj_outside
         COMMAND "${CMAKE_BINARY_DIR}/tools/line_proj" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/latlon4x2.vtk"
                 "-v" "edge_integrated_velocity"
                 "-N" "1"
                 "-p" "(-1., 22.5), (-22.5, 20.9410204722438422)"
                 "-P" "360."
                 "-verbose")
set_tests_properties(line_proj_outside PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total integrated length for segment")

add_test(NAME flux_line_2pts
         COMMAND "${CMAKE_BINARY_DIR}/tools/flux" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/cs_64.vtk"
                 "-v" "edge_integrated_velocity"
                 "-nline" "2"
                 "-xline" "0.0 + 360.*t"
                 "-yline" "-88. + 176.*t")
set_tests_properties(flux_line_2pts PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total t")

add_test(NAME flux_line_3pts
         COMMAND "${CMAKE_BINARY_DIR}/tools/flux" 
                 "-i" "${CMAKE_SOURCE_DIR}/data/cs_64.vtk"
                 "-v" "edge_integrated_velocity"
                 "-nline" "3"
                 "-xline" "0.0 + 360.*t"
                 "-yline" "-88. + 176.*t")
set_tests_properties(flux_line_3pts PROPERTIES
                     FAIL_REGULAR_EXPRESSION "Warning: total t")


add_test(NAME regrid_axis
         COMMAND "testRegridAxis")
set_tests_properties(regrid_axis PROPERTIES
                     PASS_REGULAR_EXPRESSION "Success")