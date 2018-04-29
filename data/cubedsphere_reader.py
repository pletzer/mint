import netCDF4
import numpy
import vtk

class CubedsphereReader:

    def __init__(self, filename):
        """
        Constructor
        @param filename UGRID file 
        """
        
        # read UGRID file
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

        ncells = connectivity.shape[0]

        # construct the unstructured grid as a collection of 
        # 2D cells. Each cell has its own cooordinates. Make
        # sure each cell's area is positive in lat-lon space
        # build unstructured grid
        self.vtk = {
            'points': vtk.vtkPoints(),
            'grid': vtk.vtkUnstructuredGrid()
        }

        self.vtk['points'].SetNumberOfPoints(4 * ncells)

        self.vtk['grid'].Allocate(1, 1)
        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(4)
        for icell in range(ncells):

            i00, i10, i11, i01 = connectivity[icell, :] - 1 # zero based indexing

            lon00, lat00 = lons[i00], lats[i00]
            lon10, lat10 = lons[i10], lats[i10]
            lon11, lat11 = lons[i11], lats[i11]
            lon01, lat01 = lons[i01], lats[i01]

            area013 = 0.5*( (lon10 - lon00)*(lat01 - lat00) - (lat10 - lat00)*(lon01 - lon00) )
            area231 = 0.5*( (lon01 - lon11)*(lat10 - lat11) - (lat01 - lat11)*(lon10 - lon11) )

            if area013 < 0.:
                print '--- icell = ', icell, ' area013 = ', area013, ' lon-lats = ', (lon00, lat00), (lon10, lat10), (lon01, lat01)
            if area231 < 0.:
                print '--- icell = ', icell, ' area231 = ', area231, ' lon-lats = ', (lon11, lat11), (lon01, lat01), (lon10, lat10)

            k = 4*icell

            self.vtk['points'].InsertPoint(k + 0, lon00, lat00, 0.)
            self.vtk['points'].InsertPoint(k + 1, lon10, lat10, 0.)
            self.vtk['points'].InsertPoint(k + 2, lon11, lat11, 0.)
            self.vtk['points'].InsertPoint(k + 3, lon01, lat01, 0.)

            ptIds.SetId(0, k + 0)
            ptIds.SetId(1, k + 1)
            ptIds.SetId(2, k + 2)
            ptIds.SetId(3, k + 3)
            self.vtk['grid'].InsertNextCell(vtk.VTK_QUAD, ptIds)

        self.vtk['grid'].SetPoints(self.vtk['points'])

        
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