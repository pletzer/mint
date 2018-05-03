import netCDF4
import numpy
import vtk

class UgridReader:

    TWOPI = 2. * numpy.pi

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

        pointArray = numpy.zeros((4 * ncells, 3))

        pointData = vtk.vtkDoubleArray()
        pointData.SetNumberOfComponents(3)
        pointData.SetNumberOfTuples(4 * ncells)
        pointData.SetVoidArray(pointArray, 4 * ncells * 3, 1)

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(4 * ncells)
        points.SetData(pointData)

        grid = vtk.vtkUnstructuredGrid()
        grid.Allocate(ncells, 1)
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

            if area013 < 0. or area231 < 0.:
                # this cell straddles the dateline
                # base longitude is lon00, add/remove 2*pi to reduce the cell deltas
                index10 = numpy.argmin([abs(lon10 - self.TWOPI - lon00), abs(lon10 - lon00), abs(lon10 + self.TWOPI - lon00)])
                index11 = numpy.argmin([abs(lon11 - self.TWOPI - lon00), abs(lon11 - lon00), abs(lon11 + self.TWOPI - lon00)])
                index01 = numpy.argmin([abs(lon01 - self.TWOPI - lon00), abs(lon01 - lon00), abs(lon01 + self.TWOPI - lon00)])

                lon10 += (index10 - 1) * self.TWOPI
                lon11 += (index11 - 1) * self.TWOPI
                lon01 += (index01 - 1) * self.TWOPI

            k0 = 4*icell
            k1, k2, k3 = k0 + 1, k0 + 2, k0 + 3 

            # storing coords as lon, lat, 0
            pointArray[k0, :] = lon00, lat00, 0.
            pointArray[k1, :] = lon10, lat10, 0.
            pointArray[k2, :] = lon11, lat11, 0.
            pointArray[k3, :] = lon01, lat01, 0.

            ptIds.SetId(0, k0)
            ptIds.SetId(1, k1)
            ptIds.SetId(2, k2)
            ptIds.SetId(3, k3)
            grid.InsertNextCell(vtk.VTK_QUAD, ptIds)


        grid.SetPoints(points)

        # add a cell locator
        loc = vtk.vtkCellLocator()
        loc.SetDataSet(grid)
        loc.BuildLocator()

        # store
        self.vtk = {
            'pointArray': pointArray,
            'pointData': pointData,
            'points': points,
            'grid': grid,
            'locator': loc,
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


###############################################################################

def main():
    import argparse
    from math import pi

    parser = argparse.ArgumentParser(description='Read ugrid file')
    parser.add_argument('-i', dest='input', default='', help='Specify input file')
    parser.add_argument('-V', dest='vtk_file', default='', help='Save grid in VTK file')
   
    args = parser.parse_args()

    ur = UgridReader(filename=args.input)

    if args.vtk_file:
        ur.saveToVtkFile(args.vtk_file)

if __name__ == '__main__':
    main()