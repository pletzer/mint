include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)

add_executable(latlon latlon.cxx)
target_link_libraries(latlon 
                      mint 
                      ${NETCDF_LIBRARY} 
                      ${HDF5_HL_LIBRARIES}
                      ${HDF5_LIBRARIES} 
                      ${CURL_LIBRARIES})


add_executable(regrid_edges regrid_edges.cxx)
target_link_libraries(regrid_edges
                      mint 
                      ${VTK_LIBRARIES})
