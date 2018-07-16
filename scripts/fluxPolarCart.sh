# compute for singular, radial field

# generate the data 
python ../tools/polarCart.py -nx 21 -ny 21 -o polarCart.vtk

# generate a family of contours
nt=16
a=0.3
n="10"
x0init="-0.4"
x0finl="0.05"
y0init="-0.45"
y0finl="0.02"
dx=$(python -c "print ($x0finl - $x0init)/float($n)")
dy=$(python -c "print ($y0finl - $y0init)/float($n)")
x0="$x0init"
y0="$y0init"
x0s="["
y0s="["
fluxes="["
while [ $(python -c "print ($x0 >= $x0init) and ($x0 < $x0finl)") == "True" ]; do
	x0=$(python -c "print $x0 + $dx")
	y0=$(python -c "print $y0 + $dy")
	xy=$(python generate_circle_lines.py -a $a -n $nt -x0 $x0 -y0 $y0)
	x0s="$x0s $x0,"
	y0s="$y0s $y0,"
	# compute flux
	echo "xy = $xy"
	../build/tools/line_proj -i ../tools/polarCart.vtk -p "$xy" >& res.txt
	flux=$(cat res.txt | grep total | awk '{print $4}')
	echo "flux = $flux"
	fluxes="$fluxes $flux,"
done
x0s="$x0s ]"
y0s="$y0s ]"
fluxes="$fluxes ]"

echo "a = $a nt = $nt"
echo "x0s = $x0s"
echo "y0s = $y0s"
echo "fluxes = $fluxes"

