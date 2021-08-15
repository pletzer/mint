from matplotlib import pylab
from math import sqrt
from matplotlib import rc
#rc('text', usetex=True)

ncells = [4050, 16200, 64800, 259200, 1036800, 4147200, 16588800, 66355200 ]
loopVertsErrors = [0.432506323032, 0.213229090732, 0.109496002478, 0.0546506820968, 0.0258635444994, 0.0137091999851, 0.00668680899773, 0.00435628670767 ]
loopEdgesErrors = [5.25621e-16, 7.26849e-16, 1.98279e-15, 1.97758e-15, 3.93782e-15, 3.71491e-15, 5.238e-15, 2.39574e-14 ]


pylab.loglog(ncells, loopVertsErrors, 'c--')
pylab.loglog(ncells, loopEdgesErrors, 'm-')
pylab.loglog([4000, 1e8], [1./sqrt(4000), 1./sqrt(1e8)], 'k-.')
pylab.legend(['bilinear', 'mimetic', '$N^{-1/2}$'], loc='center right')
pylab.loglog(ncells, loopVertsErrors, 'ko', markerfacecolor='None')
pylab.loglog(ncells, loopEdgesErrors, 'ko', markerfacecolor='None')
pylab.ylabel('Max loop integral error')
pylab.xlabel('Number of source grid cells')
pylab.title('Uniform to cubed sphere $6x64^2$ regridding')
pylab.show()
