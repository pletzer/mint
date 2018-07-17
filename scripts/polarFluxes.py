from matplotlib import pylab

a = 0.3; nt = 16
x0s = [ -0.355, -0.31, -0.265, -0.22, -0.175, -0.13, -0.085, -0.04, 0.005, 0.05, ]
y0s = [ -0.403, -0.356, -0.309, -0.262, -0.215, -0.168, -0.121, -0.074, -0.027, 0.02, ]
fluxes = [ 1.5662e-12, 1.37218e-12, 0.00857159, 0.189808, 0.577754, 0.936307, 1, 1, 1, 1, ]
fluxesTrapezoidal = [ 0.000035393937, 0.000237118514, 0.001934941310, 0.023454251192, 0.898960990366, 1.015476246881, 1.013067987470, 1.013052368490, 1.013052368339, 1.013052368337, ]

x = [i for i in range(len(x0s))]
pylab.plot(x, fluxes, 'm-', linewidth=2)
pylab.plot(x, fluxesTrapezoidal, '--', color='brown', linewidth=2)
pylab.legend(['mimetic', 'trapezoidal'])
pylab.plot(x, fluxes, 'ko', markerfacecolor='None')
pylab.plot(x, fluxesTrapezoidal, 'ko', markerfacecolor='None')
pylab.plot([0, 9], [1, 1], 'k-')
pylab.plot([0, 9], [0, 0], 'k-')
pylab.xlabel('Contour E index')
pylab.ylabel('flux')
pylab.title('Flux integral vs circle centre')
pylab.show()
