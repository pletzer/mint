include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)

link_directories(${NETCDF_DIR}/lib)

add_executable(regrid_edges regrid_edges.cxx)
target_link_libraries(regrid_edges
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(regrid_edges3d regrid_edges3d.cxx)
target_link_libraries(regrid_edges3d
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(line_proj line_proj.cxx)
target_link_libraries(line_proj
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(line_proj3d line_proj3d.cxx)
target_link_libraries(line_proj3d
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(flux flux.cxx)
target_link_libraries(flux
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)
