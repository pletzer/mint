include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${NETCDF_DIR}/include
)

link_directories(${NETCDF_DIR}/lib)

add_executable(latlon latlon.cxx)
target_link_libraries(latlon 
                      mint
                      netcdf
                      hdf5_hl
                      ${HDF5_LIBRARIES}
                      ${CURL_LIBRARIES}
)


add_executable(regrid_edges regrid_edges.cxx)
target_link_libraries(regrid_edges
                      mint
                      ${VTK_LIBRARIES}
                      netcdf
                      hdf5_hl
                      ${HDF5_LIBRARIES}
                      ${CURL_LIBRARIES}
)

add_executable(line_proj line_proj.cxx)
target_link_libraries(line_proj
                      mint
                      ${VTK_LIBRARIES}
                      netcdf
                      hdf5_hl
                      ${HDF5_LIBRARIES}
                      ${CURL_LIBRARIES}
)
