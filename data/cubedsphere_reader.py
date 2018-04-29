import netCDF4
import numpy
import vtk

class CubedsphereReader:

    def __init__(self, filename):
        
        nc = netCDF4.Dataset(filename, 'r')

        self.longitude = []
        self.latitude = []
        self.cell_connectivity = []

        for varname in nc.variables:
            var = nc.variables[varname]
            if hasattr(var, 'cf_role') and var.cf_role == 'face_node_connectivity':
                self.cell_connectivity = var[:]
            elif hasattr(var, 'standard_name'):
                if var.standard_name == 'longitude':
                    self.longitude = var[:]
                elif var.standard_name == 'latitude':
                    self.latitude = var[:]



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