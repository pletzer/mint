include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)

configure_file(advect_shock.sh.in advect_shock.sh @ONLY)
configure_file(advect_latlon.sh.in advect_latlon.sh @ONLY)
configure_file(advect_cubedsphere.sh.in advect_cubedsphere.sh @ONLY)
configure_file(advect_cubedsphere_periodic.sh.in advect_cubedsphere_periodic.sh @ONLY)
configure_file(advect_cubedsphere2.sh.in advect_cubedsphere2.sh @ONLY)
configure_file(advect_cubedsphere_64.sh.in advect_cubedsphere_64.sh @ONLY)

link_directories(${NETCDF_DIR}/lib)

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
                      ${NETCDF_LIBRARIES}
)

add_executable(flux flux.cxx)
target_link_libraries(flux
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)

add_executable(polySegmentIter polySegmentIter.cxx)
target_link_libraries(polySegmentIter
                      mint
                      ${VTK_LIBRARIES}
                      ${NETCDF_LIBRARIES}
)
