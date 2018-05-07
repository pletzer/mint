import netCDF4
import numpy
import vtk

class ReaderBase(object):

    TWOPI = 2. * numpy.pi

    def __init__(self):
        """
        Constructor
        No args
        """

        self.vtk = {
            'pointArray': [],
            'pointData': vtk.vtkDoubleArray(),
            'points': vtk.vtkPoints(),
            'grid': vtk.vtkUnstructuredGrid(),
            'locator': vtk.vtkCellLocator(),
        }


    def getLonLatPoints(self):
        """
        Get the longitudes and latitudes in radian at the cell vertices
        @return array
        """
        return self.vtk['pointArray'][:, :2]

        
    def saveToVtkFile(self, filename):
        """
        Save the grid to a VTK file
        @param filename VTK file
        """
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.vtk['grid'])
        writer.Update()



    def getUnstructuredGrid(self):
        """
        Get the unstructured grid
        @return vtkUnstructuredGrid instance
        """
        return self.vtk['grid']


    def getUnstructuredGridCellLocator(self):
        """
        Get the unstructured grid cell locator
        @return vtkCellLocator instance
        """    	
    	return self.vtk['locator']

