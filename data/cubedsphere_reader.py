import netCDF4
import numpy
import vtk

class CubedsphereReader:

    def __init__(self, filename):
        """
        Constructor
        @param filename UGRID file 
        """
        
        nc = netCDF4.Dataset(filename, 'r')

        self.cell_connectivity = []

        lats, lons = None, None
        connectivity = None
        for varname in nc.variables:
            var = nc.variables[varname]
            if hasattr(var, 'cf_role') and var.cf_role == 'face_node_connectivity':
                connectivity = var[:]
            elif hasattr(var, 'standard_name'):
                if var.standard_name == 'longitude':
                    lons = var[:]
                elif var.standard_name == 'latitude':
                    lats = var[:]

        nverts = len(lons)
        self.xyz = numpy.zeros((nverts, 3), numpy.float64)
        self.xyz[:, 0] = lons
        self.xyz[:, 1] = lats

        ncells = connectivity.shape[0]

        self.connectivity = numpy.zeros((ncells, 5), numpy.int64)
        self.connectivity[:, 0] = 4 # points per cell
        self.connectivity[:, 1:] = connectivity - 1 # zero based

        # build unstructured grid
        self.vtk = {
            'point_array': vtk.vtkDoubleArray(),
            'points': vtk.vtkPoints(),
            'cell_array': vtk.vtkIdTypeArray(),
            'cells': vtk.vtkCellArray(),
            'grid': vtk.vtkUnstructuredGrid()
        }

        self.vtk['point_array'].SetNumberOfComponents(3)
        self.vtk['point_array'].SetNumberOfTuples(nverts)
        self.vtk['point_array'].SetVoidArray(self.xyz, 3*nverts, 1)


        self.vtk['points'].SetNumberOfPoints(nverts)
        self.vtk['points'].SetData(self.vtk['point_array'])

        self.vtk['cell_array'].SetNumberOfComponents(5)
        self.vtk['cell_array'].SetNumberOfTuples(ncells)
        self.vtk['cell_array'].SetVoidArray(self.connectivity, 5*ncells, 1)

        self.vtk['cells'].SetNumberOfCells(ncells)
        self.vtk['cells'].SetCells(ncells, self.vtk['cell_array'])

        self.vtk['grid'].SetPoints(self.vtk['points'])

        print self.vtk['cells']
        self.vtk['grid'].Allocate(ncells, 1)
        self.vtk['grid'].SetCells(vtk.VTK_QUAD, self.vtk['cells'])


        
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

###############################################################################

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Read cubedsphere file')
    parser.add_argument('-i', dest='input', default='', help='Specify input file')
    parser.add_argument('-V', dest='vtk_file', default='', help='Save grid in VTK file')
    args = parser.parse_args()

    csr = CubedsphereReader(filename=args.input)
    if args.vtk_file:
        csr.saveToVtkFile(args.vtk_file)

if __name__ == '__main__':
    main()