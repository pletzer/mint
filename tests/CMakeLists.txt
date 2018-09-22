include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_Fortran_MODULE_DIRECTORY}
)

configure_file(testGrid.cpp testGrid.cxx)

add_executable(testLineLineIntersector testLineLineIntersector.cxx)
target_link_libraries(testLineLineIntersector
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
                     )

add_executable(testCellLocator testCellLocator.cxx)
target_link_libraries(testCellLocator
                      mint
                      ${VTK_LIBRARIES}
                      )

add_executable(testPrintAddressF testPrintAddressF.f90)
set_target_properties(testPrintAddressF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testPrintAddressF
                      mint
                      ${VTK_LIBRARIES}
                      )

add_executable(testCellLocatorF testCellLocatorF.f90)
set_target_properties(testCellLocatorF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testCellLocatorF
                      mint
                      ${VTK_LIBRARIES}
                      )

add_executable(testCellLocatorFromFileF testCellLocatorFromFileF.f90)
set_target_properties(testCellLocatorFromFileF PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(testCellLocatorFromFileF
                      mint
                      ${VTK_LIBRARIES}
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

add_test(NAME lineLineIntersector COMMAND testLineLineIntersector)
add_test(NAME polysegmentIter COMMAND testPolysegmentIter)
add_test(NAME grid COMMAND testGrid)
add_test(NAME cellLocator COMMAND testCellLocator)
add_test(NAME cellLocatorF COMMAND testCellLocatorF)
add_test(NAME cellLocatorFromFile_cs_64 COMMAND testCellLocatorFromFileF "-v" "-i" "${CMAKE_SOURCE_DIR}/data/cs_64.vtk" "-n" "10" "-o" "out.vtk")
add_test(NAME cellLocatorFromFile_lfric_24576cells COMMAND testCellLocatorFromFileF "-i" "${CMAKE_SOURCE_DIR}/data/lfric_grid.vtk" "-n" "256")
add_test(NAME cellLocatorFromFile_lfric_24576cells_check COMMAND testCellLocatorFromFileF "-i" "${CMAKE_SOURCE_DIR}/data/lfric_grid.vtk" "-n" "256" "-c")
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
add_test(NAME testPrintAddressF
         COMMAND testPrintAddressF)
add_test(NAME c2fString
         COMMAND testC2fString)

