include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
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

add_executable(testCellLocatorF testCellLocatorF.f90)
target_link_libraries(testCellLocatorF
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



add_test(NAME lineLineIntersector COMMAND testLineLineIntersector)
add_test(NAME polysegmentIter COMMAND testPolysegmentIter)
add_test(NAME grid COMMAND testGrid)
add_test(NAME cellLocator COMMAND testCellLocator)
add_test(NAME cellLocatorF COMMAND testCellLocatorF)
add_test(NAME exprParser COMMAND testExprParser)
add_test(NAME size COMMAND testSize)
