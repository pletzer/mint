cs_file_dir="$HOME/lfric-8-Jun-2018/r14523_cell_locator/mesh_tools/example/"

stream_function="sin(6*pi*x/180. - 0.3)**2 + cos(2*pi*y/180.)**2"

# d psi/d theta
u_function="-2*cos(2*pi*y/180.)*sin(2*pi*y/180.)*2*pi/180."

# -d psi/d lambda
v_function="-2*sin(6*pi*x/180. - 0.3)*cos(6*pi*x/180. - 0.3)*6*pi/180."

for n in 4 8; do # 16 32 64 128 256 512 1024 2048; do
	echo $n $stream_function $u_function $v_function
	python ../tools/ugrid_reader.py -i "${cs_file_dir}/cs_${n}.nc" \
         -stream "$stream_function" \
         -u "$u_function" -v "$v_function" \
         -V "cs_${n}.vtk"
done
