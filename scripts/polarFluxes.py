from matplotlib import pylab
import numpy

a = 0.3; nt = 16
x0s = [ -0.355, -0.31, -0.265, -0.22, -0.175, -0.13, -0.085, -0.04, 0.005, 0.05, ]
y0s = [ -0.403, -0.356, -0.309, -0.262, -0.215, -0.168, -0.121, -0.074, -0.027, 0.02, ]
fluxes = [ 1.5662e-12, 1.37218e-12, 0.00857159, 0.189808, 0.577754, 0.936307, 1, 1, 1, 1, ]
fluxesTrapezoidal = [ 0.000035393937, 0.000237118514, 0.001934941310, 0.023454251192, 0.898960990366, 1.015476246881, 1.013067987470, 1.013052368490, 1.013052368339, 1.013052368337, ]
fluxesTrapezoidalBilinear = [ -0.000046434194, -0.000310618690, -0.002464175687, -0.006981368757, 0.888861126223, 0.972763572732, 0.974484343570, 0.974495358298, 0.974495358404, 0.974495358406, ]

x = [i for i in range(len(x0s))]
pylab.plot(x, fluxes, 'm-', linewidth=4)
#pylab.plot(x, fluxesTrapezoidal, '--', color='peru', linewidth=2)
pylab.plot(x, fluxesTrapezoidalBilinear, ':', color='peru', linewidth=2)
pylab.legend(['mimetic', 'trapezoidal-exact', 'trapezoidal-bilinear'], fontsize='large')
pylab.legend(['mimetic', 'trapezoidal-bilinear'], fontsize='large')
pylab.plot(x, fluxes, 'ko', markerfacecolor='None')
#pylab.plot(x, fluxesTrapezoidal, 'ko', markerfacecolor='None')
pylab.plot(x, fluxesTrapezoidalBilinear, 'ko', markerfacecolor='None')
pylab.plot([0, 9], [1, 1], '-', color='grey')
pylab.plot([0, 9], [0, 0], '-', color='grey')
pylab.xticks( numpy.linspace(0, 9, 10) )
pylab.xlabel('Contour E index')
pylab.ylabel('Flux')
pylab.title('Flux integral vs contour')
pylab.show()
