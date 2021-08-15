from matplotlib.pylab import loglog, show, xlabel, ylabel, title
import sys

num_cells = eval(sys.argv[1])
relerrors = eval(sys.argv[2])

loglog(num_cells, relerrors, 'm-', num_cells, relerrors, 'ko', markerfacecolor='None')

alpha = -1. # convergence rate in number of cells
n0, n1 = num_cells[0], num_cells[-1]
e0 = 1.0
c = n0**(-alpha) * e0
e1 = c * n1**alpha
loglog([n0, n1], [e0, e1], 'k--')

xlabel('Number of horizontal grid cells')
ylabel('Relative error')
title('Flux error on cubed sphere')
show()

