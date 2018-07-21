from matplotlib import pylab
from math import sqrt

ncells = [4050, 16200, 64800, 259200, 1036800, 4147200, 16588800, 66355200 ]

loopVertsErrors = [0.432506323032, 0.213229090732, 0.109496002478, 0.0546506820968, 0.0258635444994, 0.0137091999851, 0.00668680899773, 0.00435628670767 ]
loopEdgesErrors = [5.25621e-16, 7.26849e-16, 2.41231e-14, 2.35502e-14, 2.32577e-14, 3.83547e-14, 3.83174e-14,  float('nan')]

pylab.loglog(ncells, loopVertsErrors, 'c--')
pylab.loglog(ncells, loopEdgesErrors, 'm-')
pylab.legend(['bilinear', 'mimetic'], loc='center right')
pylab.loglog(ncells, loopVertsErrors, 'ko', markerfacecolor='None')
pylab.loglog(ncells, loopEdgesErrors, 'ks', markerfacecolor='None')
pylab.loglog([4000, 1e8], [1./sqrt(4000), 1./sqrt(1e8)], 'k-.')
pylab.ylabel('Max loop integral error')
pylab.xlabel('Number of source grid cells')
pylab.title('Uniform to cubed sphere $6x64^2$ regridding')
pylab.show()
