 python latlon.py -nlon 100 -nlat 60 -o um100x60.nc
 python ugrid_reader.py -i ../data/cubedsphere10.nc -V cs.vtk
 python latlon_reader.py -i um100x60.nc -stream "x*cos(y)" -p 30 -V um100x60.vtk
 python regrid_edges.py -s um100x60.vtk -v "edge_integrated_velocity" -d cs.vtk -o regrid.vtk

