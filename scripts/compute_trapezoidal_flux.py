import vtk
import argparse
import numpy

parser = argparse.ArgumentParser(description='Write line point data to file')
parser.add_argument('-p', type=str, default="(0., 0.),(1., 0.)", 
                    help='Interlaced xy points')
parser.add_argument('-o', type=str, default="line.vtk", 
                    help='Output file')
args = parser.parse_args()

xy = numpy.array(eval('[' + args.p + ']'))
npts = xy.shape[0]
ncells = npts - 1
twopi = 2.*numpy.pi

def vectorField(xy):
	x, y = xy
	r2 = x**2 + y**2
	return numpy.array([x, y, 0.])/(twopi*r2)

def integratedFlux(xy0, xy1):
	xyMid = 0.5*(xy0 + xy1)
	vec = vectorField(xyMid)
	ds = numpy.array([xy1[1] - xy0[1], -(xy1[0] - xy0[0]), 0.0])
	return numpy.dot(vec, ds)

flux = 0.0
for i in range(ncells):
	xy0 = xy[i + 0, :]
	xy1 = xy[i + 1, :]
	flux += integratedFlux(xy0, xy1)

print 'flux = {:16.12f}'.format(flux)