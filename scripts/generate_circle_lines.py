import numpy
import argparse

parser = argparse.ArgumentParser(description='Generate line point data')
parser.add_argument('-n', type=int, default=32, 
                    help='Number of points')
parser.add_argument('-a', type=float, default=0.25, 
                    help='Radius')
parser.add_argument('-x0', type=float, default=0., 
                    help='Radius centre x coordinate')
parser.add_argument('-y0', type=float, default=0., 
                    help='Radius centre y coordinate')
args = parser.parse_args()

ts = numpy.linspace(0., 2*numpy.pi, args.n + 1)
xs = args.x0 + args.a*numpy.cos(ts)
ys = args.y0 + args.a*numpy.sin(ts)

res = ''
for i in range(args.n + 1):
	res += '({:15.12f}, {:15.12f}),'.format(xs[i], ys[i])
print res
