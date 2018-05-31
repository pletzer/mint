include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)

add_executable(testLineLineIntersector testLineLineIntersector.cxx)
target_link_libraries(testLineLineIntersector
	                  mint 
	                  ${NETCDF_LIBRARY} 
	                  ${HDF5_HL_LIBRARIES}
	                  ${HDF5_LIBRARIES} 
	                  ${CURL_LIBRARIES})

