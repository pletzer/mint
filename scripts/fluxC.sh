# compute flux

# function determining the path
nline=2
xlinefunc="360*t"
ylinefunc="0.*t -35.0"

# write to VTK file
python writeLineToVTK.py "$nline" "$xlinefunc" "$ylinefunc" lineC.vtk

# stream function
streamfunc="sin(2*pi*x/180. - 0.3)**2 + cos(2*pi*y/180. - 0.2)**2"

x0expr=$(echo $xlinefunc | perl -ne "s/t/0./g;print;")
y0expr=$(echo $ylinefunc | perl -ne "s/t/0./g;print;")
x1expr=$(echo $xlinefunc | perl -ne "s/t/1./g;print;")
y1expr=$(echo $ylinefunc | perl -ne "s/t/1./g;print;")

x0=$(python -c "from math import *; print '{:3.18f}'.format($x0expr)")
y0=$(python -c "from math import *; print '{:3.18f}'.format($y0expr)")
x1=$(python -c "from math import *; print '{:3.18f}'.format($x1expr)")
y1=$(python -c "from math import *; print '{:3.18f}'.format($y1expr)")

# exact flux
s0=$(python -c "from math import *; x = $x0; y = $y0; print '{:3.18f}'.format($streamfunc)")
s1=$(python -c "from math import *; x = $x1; y = $y1; print '{:3.18f}'.format($streamfunc)")
exact_flux=$(python -c "print '{:3.18f}'.format($s1 - $s0)")
echo "Exact flux = $exact_flux"


fluxes="["
num_cells="["
errors="["
relerrors="["
for n in 4 8 16 32 64 128 256 512 1024 2048; do
	echo "n = $n"

	rm -rf rs.txt
	../build/tools/flux -i cs_${n}.vtk \
	             -nline "$nline" -xline "$xlinefunc" -yline "$ylinefunc" \
	              >& res.txt

	ncells=$(python get_num_cells.py "res.txt")
	flux=$(python get_flux.py "res.txt")
	err=$(python -c "print '{:3.18f}'.format(abs($flux - $exact_flux))")
	relerr=$(python -c "print '{:3.18f}'.format($err/abs($exact_flux))")

	echo "Flux error: $err"

	num_cells="$num_cells $ncells,"
	fluxes="$fluxes $flux,"
	errors="$errors $err,"
	relerrors="$relerrors $relerr,"
done
num_cells="$num_cells ]"
fluxes="$fluxes ]"
errors="$errors ]"
relerrors="$relerrors ]"

# 
echo "num_cells  = $num_cells"
echo "fluxes     = $fluxes"
echo "errors     = $errors"
echo "rel_errors = $relerrors"

# plot
python plot.py "$num_cells" "$relerrors"
