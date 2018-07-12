import numpy
import math
import argparse
import vtk
from math import cos, sin, pi, log, exp

"""
A class to generate a polar grid
"""

class Polar:

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
            'velocity': vtk.vtkDoubleArray(),
        }



    def setNumberOfRhoCells(self, n):
        """
        Set the number of radial cells
        @param n number of cells
        """
        self.numRho0 = n

    def setNumberOfTheCells(self, n):
        """
        Set the number of poloidal cells
        @param n number of cells
        """
        self.numThe0 = n

    def build(self, potentialExpr, radius=1.0):
        """
        Build the object. Call this after setNumberOfRhoCells and setNumberOfTheCells
        """

        # deltas
        self.dRho = radius/float(self.numRho0)
        self.dThe = 2. * math.pi/float(self.numThe0)

        # axes
        self.rhos = numpy.linspace(0.0, radius, self.numRho0 + 1)
        self.thes = numpy.linspace(0.0, 2*math.pi, self.numThe0 + 1)

        # construct the unstructured grid as a collection of 
        # 2D cells. Each cell has its own cooordinates. Make
        # sure each cell's area is positive in lat-lon space
        # build unstructured grid


        ncells = self.numRho0 * self.numThe0
        pointArray = numpy.zeros((4 * ncells, 3))
        self.vtk['pointArray'] = pointArray

        edgeData = self.vtk['velocity']
        edgeData.SetNumberOfComponents(4)
        edgeData.SetNumberOfTuples(ncells)
        edgeData.SetName('edgeData')

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

        icell = 0
        for i in range(self.numRho0):
            rho0 = self.rhos[i]
            rho1 = self.rhos[i + 1]

            for j in range(self.numThe0):
                the0 = self.thes[j]
                the1 = self.thes[j + 1]

                x00, y00 = rho0*math.cos(the0), rho0*math.sin(the0)
                x10, y10 = rho1*math.cos(the0), rho1*math.sin(the0)
                x11, y11 = rho1*math.cos(the1), rho1*math.sin(the1)
                x01, y01 = rho0*math.cos(the1), rho0*math.sin(the1)

                x, y = x00, y00
                psi00 = eval(potentialExpr)
                x, y = x10, y10
                psi10 = eval(potentialExpr)
                x, y = x11, y11
                psi11 = eval(potentialExpr)
                x, y = x01, y01
                psi01 = eval(potentialExpr)

                edgeData.SetTuple(icell, [psi10 - psi00, psi11 - psi10, psi01 - psi11, psi00 - psi01])

                k0 = 4*icell
                k1, k2, k3 = k0 + 1, k0 + 2, k0 + 3 

                pointArray[k0, :] = x00, y00, 0.
                pointArray[k1, :] = x10, y10, 0.
                pointArray[k2, :] = x11, y11, 0.
                pointArray[k3, :] = x01, y01, 0.

                ptIds.SetId(0, k0)
                ptIds.SetId(1, k1)
                ptIds.SetId(2, k2)
                ptIds.SetId(3, k3)
                grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

                icell += 1

        grid.SetPoints(points)
        grid.GetCellData().AddArray(edgeData)


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
    parser.add_argument('-nrho', dest='nrho', type=int, default=10, help='Number of rho cells')
    parser.add_argument('-nthe', dest='nthe', type=int, default=15, help='Number of the cells')
    parser.add_argument('-o', dest='output', type=str, default='polar.vtk', help='Output file')
    parser.add_argument('-stream', dest='streamFunc', default='0.5*log(x**2 + y**2 + 1.e-14)/(2.*pi)', help='')


    args = parser.parse_args()
    numRho, numThe = args.nrho, args.nthe
    outputFile = args.output

    pl = Polar()
    pl.setNumberOfRhoCells(numRho)
    pl.setNumberOfTheCells(numThe)
    pl.build(args.streamFunc)
    pl.saveToVtkFile(outputFile)

if __name__ == '__main__':
    main()

