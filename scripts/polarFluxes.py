from matplotlib import pylab

a = 0.3; nt = 16
x0s = [ -0.355, -0.31, -0.265, -0.22, -0.175, -0.13, -0.085, -0.04, 0.005, 0.05, ]
y0s = [ -0.403, -0.356, -0.309, -0.262, -0.215, -0.168, -0.121, -0.074, -0.027, 0.02, ]
fluxes = [ 1.5662e-12, 1.37218e-12, 0.00857159, 0.189808, 0.577754, 0.936307, 1, 1, 1, 1, ]

pylab.plot(x0s, fluxes, 'm-')
pylab.plot(x0s, fluxes, 'ko', markerfacecolor='None')
pylab.xlabel('circle centre x coordinate')
pylab.ylabel('flux')
pylab.title('Flux integral vs circle centre')
pylab.show()
