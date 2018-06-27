import sys
from numpy import *

nline = eval(sys.argv[1])
t = linspace(0., 1., nline)
x = eval(sys.argv[2])
y = eval(sys.argv[3])

f = open(sys.argv[4], 'w')
f.write('# vtk DataFile Version 4.2\n')
f.write('vtk output\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('POINTS {} double\n'.format(nline))
for i in range(nline):
	f.write('{} {} 0\n'.format(x[i], y[i]))
ncells = nline - 1
f.write('CELLS {} {}\n'.format(ncells, ncells*3))
for i in range(ncells):
	f.write('2 {} {}\n'.format(i, i + 1))
f.write('CELL_TYPES {}\n'.format(ncells))
for i in range(ncells):
	f.write('3\n')
f.close()
