import vtk
import numpy
import argparse
from scipy.integrate import odeint
import random
import functools

"""
Convert an edge field defined cell by cell to an edge field defined on a collection of line cells
"""

LON_INDEX, LAT_INDEX = 0, 1
EPS = 1.234e-12
TOL = 1.e-36 # 1.e-10
subId = vtk.mutable(-1)
xsis = numpy.zeros((3,), numpy.float64)
etas = numpy.zeros((3,), numpy.float64)
weights = numpy.zeros((8,), numpy.float64)
cell = vtk.vtkGenericCell()

parser = argparse.ArgumentParser(description='Convert cell-by-cell edge field into an edge field')
parser.add_argument('-i', dest='inputFile', default='res.vtk', help='Specify path to VTK, cell-by-cell input file')
parser.add_argument('-o', dest='outputFile', default='trajectory.vtk', help='Specify name of VTK output file')
parser.add_argument('-v', dest='edgeFieldName', default='edge_integrated_velocity', help='Specify name of edge integrated variable')
parser.add_argument('-tf', dest='finalTime', default=100.0, type=float, help='Specify final time')
parser.add_argument('-nt', dest='numSteps', default=100, type=int, help='Specify number of time steps')
parser.add_argument('-ns', default=10, type=int, 
                    help='Number of random seed points')


args = parser.parse_args()

def saveTrajectory(sols, outputFile):
    """
    Save the trajectory to VTK file
    @param sols list of return values of odeint
    @param filename
    """

    # number of contours
    nContours = len(sols)

    # number of points for each contour
    nptsContour = [sol.shape[0] for sol in sols]

    # total number of points
    npts = functools.reduce(lambda x, y: x + y, nptsContour)

    # total number of segments
    nSegs = functools.reduce(lambda x, y: x + y, [nps - 1 for nps in nptsContour])

    # number of space dimensions
    ndims = 3

    pvals = numpy.zeros((npts, 3), numpy.float64)
    tarr = vtk.vtkDoubleArray()
    tpts = vtk.vtkPoints()
    tgrid = vtk.vtkUnstructuredGrid()

    tarr.SetNumberOfComponents(ndims)
    tarr.SetNumberOfTuples(npts)
    tpts.SetNumberOfPoints(npts)

    ptIds = vtk.vtkIdList()
    ptIds.SetNumberOfIds(2)

    tgrid.Allocate(nSegs, 1)

    # create the points and the unstructured grid that goes with it
    offset1 = 0
    offset2 = 0
    for iContour in range(nContours):

        ns = nptsContour[iContour]

        # store points
        for i in range(ns):
            pvals[i + offset1, :] = sols[iContour][i]
            pvals[i + offset1, 0] = max(0., min(360., pvals[i + offset1, 0]))
            pvals[i + offset1, 1] = max(-90., min(90., pvals[i + offset1, 1]))
        offset1 += ns

        # create new cells/segments
        for i in range(ns - 1):
            ptIds.SetId(0, i + offset2)
            ptIds.SetId(1, i + 1 + offset2)
            tgrid.InsertNextCell(vtk.VTK_LINE, ptIds)
        offset2 += ns

    # connect
    tpts.SetData(tarr)
    tgrid.SetPoints(tpts)
    tarr.SetVoidArray(pvals, npts*3, 1)

    # save
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(outputFile)
    writer.SetInputData(tgrid)
    writer.Update()


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
    if data is None:
        print('ERROR: no field called {}'.format(args.edgeFieldName))
        return numpy.zeros((3,), numpy.float64)

    # apply periodicity on longitudes when leaving domain
    if point[0] < 0.:
        point[0] += 360.
    elif point[0] > 360.:
        point[0] -= 360.

    # lat cannot exceed 90 deg
    point[1] = min(90., point[1])
    point[1] = max(-90., point[1])

    # find the cell and the param coords xis
    cellId = loc.FindCell(point, TOL, cell, xsis, weights)
    if cellId < 0:
        print('ERROR: out of domain integration, point = {}'.format(point[:2]))
        return numpy.zeros((3,), numpy.float64)

    # complement to xi
    etas[:] = 1.0 - xsis

    # get the corner points of the cell
    ptIds = cell.GetPointIds()
    corners = []
    for i in range(ptIds.GetNumberOfIds()):
        p = pts.GetPoint(ptIds.GetId(i))
        corners.append( numpy.array([p[0], p[1]]) )

    # d position / d xi^0
    dPd0 = (corners[1] - corners[0])*etas[1] + (corners[2] - corners[3])*xsis[1]
    # d position / d xi^1
    dPd1 = (corners[3] - corners[0])*etas[0] + (corners[2] - corners[1])*xsis[0]

    # Jacobian (cell area in lon-lat coords)
    jac = dPd0[0]*dPd1[1] - dPd0[1]*dPd1[0]

    edgeVals = numpy.array(data.GetTuple(cellId))

    yp = numpy.zeros((3,), numpy.float64)
    yp[0] = ( - (edgeVals[0]*etas[1] + edgeVals[2]*xsis[1])*dPd1[0] 
              + (edgeVals[3]*etas[0] + edgeVals[1]*xsis[0])*dPd0[0] )/jac
    yp[1] = ( - (edgeVals[0]*etas[1] + edgeVals[2]*xsis[1])*dPd1[1] 
              + (edgeVals[3]*etas[0] + edgeVals[1]*xsis[0])*dPd0[1] )/jac

    return yp

# read the file 
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(args.inputFile)
reader.Update()

# get the unstructured grid
ugrid = reader.GetOutput()

npts = ugrid.GetNumberOfPoints()
ncells = ugrid.GetNumberOfCells()

# create a locator
loc = vtk.vtkCellLocator()
loc.SetDataSet(ugrid)
loc.BuildLocator()

timeSteps = numpy.linspace(0.0, args.finalTime, args.numSteps + 1)

sols = []
for isol in range(args.ns):
    p0 = numpy.array([0. + 360.*random.random(), -90. + 180.*random.random(), 0.0])
    sol = odeint(tendency, p0, timeSteps, tfirst=True, args=(loc, ugrid))
    sols.append(sol)

# save the trajectory
print('saving the trajectory in {}'.format(args.outputFile))
saveTrajectory(sols, args.outputFile)
