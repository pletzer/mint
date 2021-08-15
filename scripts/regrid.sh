stream_function="sin(2*pi*x/180. - 0.3)**2 + cos(2*pi*y/180. - 0.2)**2"

ncells=""
loopVertsErrors=""
loopEdgesErrors=""
for nlat in 45 90 180 360 720 1440 2880 5760; do # 11520; do
	nlon=$(expr 2 \* $nlat)
	ntot=$(expr $nlon \* $nlat)
	echo "nlon = $nlon nlat = $nlat $(date)"
	# generate the grid in UM netcdf format
	python ../tools/latlon.py -nlon $nlon -nlat $nlat -o um.nc
	# convert to VTK format and add vector data
	python ../tools/latlon_reader.py -i um.nc -stream "$stream_function" -V um.vtk -b

	# regrid using bilinear interpolation to cubed sphere 6 * 64^2 resolution
	echo "regrid_verts...$(date)"
	python ../tools/regrid_verts.py -s um.vtk -v "edge_integrated_velocity" -d cs_64.vtk -o regrid_verts.vtk >& verts.txt

	# regrid using mimetic interpolation
	echo "regrid_edges...$(date)"
	../build-lfric/tools/regrid_edges -N 10000 -s um.vtk -v "edge_integrated_velocity" -d cs_64.vtk -o regrid_edges.vtk >& edges.txt
	# extract the max loop integral error
	max_loop_error_verts=$(cat verts.txt | awk -F '/' '{print $5}')
	max_loop_error_edges=$(cat edges.txt | perl -ne 'if(/Min/){print;}' | awk -F '/' '{print $5}')
	ncells="$ncells, $ntot"
	loopVertsErrors="$loopVertsErrors, $max_loop_error_verts"
	loopEdgesErrors="$loopEdgesErrors, $max_loop_error_edges"

	echo "ncells = $ncells"
	echo "loopVertsErrors = $loopVertsErrors"
	echo "loopEdgesErrors = $loopEdgesErrors"

done
ncells="[$ncells ]"
loopVertsErrors="[$loopVertsErrors ]"
loopEdgesErrors="[$loopEdgesErrors ]"
echo "ncells = $ncells"
echo "loopVertsErrors = $loopVertsErrors"
echo "loopEdgesErrors = $loopEdgesErrors"

