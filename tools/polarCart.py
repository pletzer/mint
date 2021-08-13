import numpy
import math
import argparse
import vtk
from math import cos, sin, pi, log, exp, atan2

"""
A class to generate a polar field on a Cartesian grid
"""

class PolarCart:

    def __init__(self):
        """
        Constructor
        no arguments
        """
        self.vtk = {
            'pointArray': [],
            'pointData': vtk.vtkDoubleArray(),
            'points': vtk.vtkPoints(),
            'grid': vtk.vtkUnstructuredGrid(),
            'edgeData': vtk.vtkDoubleArray(),
            'vectorField': vtk.vtkDoubleArray(),
        }



    def setNumberOfCells(self, nx, ny):
        """
        Set the number of radial cells
        @param nx number of cells in the x direction
        @param ny number of cells in the y direction
        """
        self.nx = nx
        self.ny = ny


    def build(self, boxsize=((-1., -1.), (1., 1.))):
        """
        Build the object. Call this after setNumberOfCells
        @param boxsize ((xmin, ymin), *xmax, ymax))
        """

        xmin, ymin = boxsize[0]
        xmax, ymax = boxsize[1]

        # deltas
        dx = (xmax - xmin)/float(self.nx)
        dy = (ymax - ymin)/float(self.ny)

        # construct the unstructured grid as a collection of 
        # 2D cells. Each cell has its own coordinates. Make
        # sure each cell's area is positive in lat-lon space
        # build unstructured grid


        ncells = self.nx * self.ny
        pointArray = numpy.zeros((4 * ncells, 3))
        self.vtk['pointArray'] = pointArray

        edgeData = self.vtk['edgeData']
        edgeData.SetNumberOfComponents(4)
        edgeData.SetNumberOfTuples(ncells)
        edgeData.SetName('edgeData')

        vectorField = self.vtk['vectorField']
        vectorField.SetNumberOfComponents(3)
        vectorField.SetNumberOfTuples(ncells)
        vectorField.SetName('vectorField')

        pointData = self.vtk['pointData']
        pointData.SetNumberOfComponents(3)
        pointData.SetNumberOfTuples(4 * ncells)
        pointData.SetVoidArray(pointArray, 4 * ncells * 3, 1)

        points = self.vtk['points']
        points.SetNumberOfPoints(4 * ncells)
        points.SetData(pointData)

        grid = self.vtk['grid']
        grid.Allocate(ncells, 1)
        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(4)

        twopi = 2 * math.pi

        icell = 0
        for i in range(self.nx):
            x0 = xmin + i*dx
            x1 = x0 + dx

            for j in range(self.ny):
                y0 = ymin + j*dy
                y1 = y0 + dy

                # difference in poloidal angle gives flux across theta edge
                # assume the radial field starts at (0, 0)
                t00 = atan2(y0, x0)
                t10 = atan2(y0, x1)
                t11 = atan2(y1, x1)
                t01 = atan2(y1, x0)

                # add/remove 2*pi 
                dts0arr = [t10 - t00 - twopi, t10 - t00, t10 - t00 + twopi]
                dt1sarr = [t11 - t10 - twopi, t11 - t10, t11 - t10 + twopi]
                dts1arr = [t11 - t01 - twopi, t11 - t01, t11 - t01 + twopi]
                dt0sarr = [t01 - t00 - twopi, t01 - t00, t01 - t00 + twopi]

                dts0 = dts0arr[numpy.argmin(numpy.fabs(dts0arr))] / twopi
                dt1s = dt1sarr[numpy.argmin(numpy.fabs(dt1sarr))] / twopi
                dts1 = dts1arr[numpy.argmin(numpy.fabs(dts1arr))] / twopi
                dt0s = dt0sarr[numpy.argmin(numpy.fabs(dt0sarr))] / twopi

                edgeData.SetTuple(icell, [dts0, dt1s, dts1, dt0s])
                vx = 0.5*(dts0 + dts1) / dx
                vy = 0.5*(dt0s + dt1s) / dy
                # grad psi cross z hat
                vectorField.SetTuple(icell, [vy, -vx, 0.])

                k0 = 4*icell
                k1, k2, k3 = k0 + 1, k0 + 2, k0 + 3 

                pointArray[k0, :] = x0, y0, 0.
                pointArray[k1, :] = x1, y0, 0.
                pointArray[k2, :] = x1, y1, 0.
                pointArray[k3, :] = x0, y1, 0.

                ptIds.SetId(0, k0)
                ptIds.SetId(1, k1)
                ptIds.SetId(2, k2)
                ptIds.SetId(3, k3)
                grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

                icell += 1

        grid.SetPoints(points)
        grid.GetCellData().AddArray(edgeData)
        grid.GetCellData().AddArray(vectorField)


    def saveToVtkFile(self, filename):
        """
        Save the grid to a VTK file
        @param filename VTK file
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.vtk['grid'])
        writer.Update()


#############################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-nx', dest='nx', type=int, default=10, help='Number of x cells')
    parser.add_argument('-ny', dest='ny', type=int, default=10, help='Number of y cells')
    parser.add_argument('-o', dest='output', type=str, default='polarCart.vtk', help='Output file')


    args = parser.parse_args()
    outputFile = args.output

    pl = PolarCart()
    pl.setNumberOfCells(args.nx, args.ny)
    pl.build()
    pl.saveToVtkFile(outputFile)

if __name__ == '__main__':
    main()

