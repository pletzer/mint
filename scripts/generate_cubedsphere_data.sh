cs_file_dir=$HOME/lfric-8-Jun-2018/r14523_cell_locator/mesh_tools/example/

for n in 4 8 16 32 64 128 256 512 1024 2048; do
	echo $n
	python ../tools/ugrid_reader.py -i ${cs_file_dir}/cs_${n}.nc \
         -stream "sin(6*pi*x/180. - 0.3)**2 + cos(2*pi*y/180.)**2" -V cs_${n}.vtk
done
