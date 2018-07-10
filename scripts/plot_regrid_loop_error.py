from matplotlib import pylab
from math import sqrt

ncells = [4050, 16200, 64800, 259200, 1036800, 4147200, 16588800, 66355200  ]
loopVertsErrors = [0.43250632328, 0.21322909097, 0.109496002498, 0.0546506820885, 0.0258448920252, 0.0137091999867, 0.00668680902321,  0.00435628670771]
loopEdgesErrors = [6.33174e-16, 3.88578e-16, 1.97238e-15, 1.98279e-15, 1.49378e-13, 1.3365e-13, 1.28068e-13,  2.70251e-13]
pylab.loglog(ncells, loopVertsErrors, 'c--')
pylab.loglog(ncells, loopEdgesErrors, 'm-')
pylab.legend(['bilinear', 'mimetic'])
pylab.loglog(ncells, loopVertsErrors, 'ko', markerfacecolor='None')
pylab.loglog(ncells, loopEdgesErrors, 'ks', markerfacecolor='None')
pylab.loglog([4000, 1e8], [1./sqrt(4000), 1./sqrt(1e8)], 'k-.')
pylab.ylabel('Max loop integral error')
pylab.xlabel('Number of source grid cells')
pylab.title('Uniform to cubed sphere $6x64^2$ regridding')
pylab.show()
