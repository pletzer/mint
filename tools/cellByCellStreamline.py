import vtk
import numpy
import argparse
from scipy.integrate import odeint

"""
Convert an edge field defined cell by cell to an edge field defined on a collection of line cells
"""

LON_INDEX, LAT_INDEX = 0, 1
EPS = 1.234e-12
TOL = 1.e-10
cell = vtk.vtkGenericCell()
xsis = numpy.zeros((3,), numpy.float64)
etas = numpy.zeros((3,), numpy.float64)
weights = numpy.zeros((8,), numpy.float64)

parser = argparse.ArgumentParser(description='Convert cell-by-cell edge field into an edge field')
parser.add_argument('-i', dest='inputFile', default='res.vtk', help='Specify path to VTK, cell-by-cell input file')
parser.add_argument('-o', dest='outputFile', default='edges.vtk', help='Specify name of VTK output file')
parser.add_argument('-v', dest='edgeFieldName', default='edge_integrated_velocity', help='Specify name of edge integrated variable')
parser.add_argument('-0', dest='initialPosition', default='180.0, 0.0', help='Specify initial condition "lon,lat"')
parser.add_argument('-tf', dest='finalTime', default=100.0, type=float, help='Specify final time')
parser.add_argument('-nt', dest='numSteps', default=100, type=int, help='Specify number of time steps')

args = parser.parse_args()



def tendency(t, point, loc, grid):
    """
    Compute the ODE tendency
    @param t time (not used)
    @param point point
    @param loc cell locator
    @param grid unstructed grid instance
    @return  velocity at the point
    """

    pts = grid.GetPoints()
    data = grid.GetCellData().GetAbstractArray(args.edgeFieldName)

    # TO DO 
    # apply periodicity on longitudes

    # find the cell and the param coords xis
    cellId = loc.FindCell(point, TOL, cell, xsis, weights)
    if cellId < 0:
        print('ERROR: out of domain integration, point = {}'.format(point[:2]))
        return numpy.zeros((3,), numpy.float64)

    # complement to xi
    etas[:] = 1.0 - xsis

    # get the corner points of the cell
    ptIds = cell.GetPointIds()
    corner = []
    for i in range(ptIds.GetNumberOfIds()):
        p = pts.GetPoint(ptIds.GetId(i))
        corner.append( numpy.array([p[0], p[1]]) )

    dPd0 = (corner[1] - corner[0])*etas[1] + (corner[2] - corner[3])*xsis[1]
    dPd1 = (corner[3] - corner[0])*etas[0] + (corner[2] - corner[1])*xsis[0]

    # Jacobian (cell area in lon-lat coords)
    jac = dPd0[0]*dPd1[1] - dPd0[1]*dPd1[0]

    edgeVals = numpy.array(data.GetTuple(cellId))

    yp = numpy.zeros((3,), numpy.float64)
    yp[0] = ( + (edgeVals[0]*etas[1] + edgeVals[2]*xsis[1])*dPd1[1] 
              - (edgeVals[3]*etas[0] + edgeVals[1]*xsis[0])*dPd0[1] )/jac
    yp[1] = ( - (edgeVals[0]*etas[1] + edgeVals[2]*xsis[1])*dPd1[0] 
              + (edgeVals[3]*etas[0] + edgeVals[1]*xsis[0])*dPd0[0] )/jac

    return yp

# read the file 
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.inputFile)
reader.Update()

# get the unstructured grid
ugrid = reader.GetOutput()

# create a locator
loc = vtk.vtkCellLocator()
loc.SetDataSet(ugrid)
loc.BuildLocator()


p0 = numpy.array(list(eval(args.initialPosition)) + [0.0])
timeSteps = numpy.linspace(0.0, args.finalTime, args.numSteps + 1)
sol = odeint(tendency, p0, timeSteps, tfirst=True, args=(loc, ugrid))

print(sol)
