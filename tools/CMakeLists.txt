include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)

link_directories(${NETCDF_DIR}/lib)

add_executable(latlon latlon.cxx)
target_link_libraries(latlon 
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(regrid_edges regrid_edges.cxx)
target_link_libraries(regrid_edges
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(line_proj line_proj.cxx)
target_link_libraries(line_proj
                      mint
                      ${VTK_LIBRARIES}
)

add_executable(stream_function_flux stream_function_flux.cxx)
target_link_libraries(stream_function_flux
                      mint
                      ${VTK_LIBRARIES}
)
