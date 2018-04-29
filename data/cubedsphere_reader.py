import netCDF4
import numpy
import vtk

class CubedsphereReader:

    def __init__(self, filename):
        
        nc = netCDF4.Dataset(filename, 'r')

        self.cell_connectivity = []

        lats, lons = None, None
        for varname in nc.variables:
            var = nc.variables[varname]
            if hasattr(var, 'cf_role') and var.cf_role == 'face_node_connectivity':
                self.cell_connectivity = var[:]
            elif hasattr(var, 'standard_name'):
                if var.standard_name == 'longitude':
                    lons = var[:]
                elif var.standard_name == 'latitude':
                    lats = var[:]

        nverts = len(lons)
        self.xyz = numpy.zeros((nverts, 3), numpy.float64)
        self.xyz[:, 0] = lats
        self.xyz[:, 1] = lons

        ncells = self.cell_connectivity.shape[0]
        # use zero-based indexing from now on
        self.cell_connectivity -= 1

        # build unstructured grid
        self.vtk = {
            'point_array': vtk.vtkDoubleArray(),
            'points': vtk.vtkPoints(),
            'cell_array': vtk.vtkIntArray(),
            'cells': vtk.vtkCellArray(),
            'grid': vtk.vtkUnstructuredGrid()
        }

        self.vtk['point_array'].SetNumberOfComponents(3)
        self.vtk['point_array'].SetNumberOfTuples(nverts)
        self.vtk['point_array'].SetVoidArray(self.xyz, 3*nverts, 1)


        self.vtk['points'].SetNumberOfPoints(nverts)
        self.vtk['points'].SetData(self.vtk['point_array'])

        self.vtk['cell_array'].SetNumberOfComponents(4)
        self.vtk['cell_array'].SetNumberOfTuples(ncells)
        self.vtk['cell_array'].SetVoidArray(self.cell_connectivity, 4*ncells, 1)

        self.vtk['grid'].SetPoints(self.vtk['points'])
        self.vtk['grid'].SetCells(vtk.VTK_QUAD, self.vtk['cells'])



        



    def getUnstructuredGrid(self):
        pass

###############################################################################

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Read cubedsphere file')
    parser.add_argument('-i', dest='input', default='', help='Specify input file')
    args = parser.parse_args()

    csr = CubedsphereReader(filename=args.input)

if __name__ == '__main__':
    main()