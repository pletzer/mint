import netCDF4
import numpy
import vtk
from reader_base import ReaderBase

class LatlonReader(ReaderBase):


    def __init__(self, filename):
        """
        Constructor
        @param filename UGRID file 
        """

        super(LatlonReader, self).__init__()
        
        # read UGRID file
        nc = netCDF4.Dataset(filename, 'r')

        # gather all the latitudes and longitudes
        lats, lons = None, None
        lats_0, lons_0 = None, None
        for varname in nc.variables:
            var = nc.variables[varname]
            if hasattr(var, 'standard_name'):
                if var.standard_name == 'longitude':
                    if varname.find('_0') >= 0:
                        lons_0 = var[:]
                    else:
                        lons = var[:]
                elif var.standard_name == 'latitude':
                    if varname.find('_0') >= 0:
                        lats_0 = var[:]
                    else:
                        lats = var[:]

        ncells_lat, ncells_lon = len(lats_0), len(lons_0)
        ncells = ncells_lat * ncells_lon

        # construct the unstructured grid as a collection of 
        # 2D cells. Each cell has its own cooordinates. Make
        # sure each cell's area is positive in lat-lon space
        # build unstructured grid

        pointArray = numpy.zeros((4 * ncells, 3))
        self.vtk['pointArray'] = pointArray

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
        for j0 in range(ncells_lat):

            j1 = j0 + 1

            for i0 in range(ncells_lon):

                i1 = (i0 + 1) % ncells_lon

                lon00, lat00 = lons[i0], lats[j0]
                lon10, lat10 = lons[i1], lats[j0]
                lon11, lat11 = lons[i1], lats[j1]
                lon01, lat01 = lons[i0], lats[j1]

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

                icell += 1


        grid.SetPoints(points)

        # build cell locator
        loc = self.vtk['locator']
        loc.SetDataSet(grid)
        loc.BuildLocator()


###############################################################################

def main():
    import argparse
    from math import pi

    parser = argparse.ArgumentParser(description='Read ugrid file')
    parser.add_argument('-i', dest='input', default='ll.nc', help='Specify input file')
    parser.add_argument('-V', dest='vtk_file', default='lonlat.vtk', help='Save grid in VTK file')
   
    args = parser.parse_args()

    lr = LatlonReader(filename=args.input)

    if args.vtk_file:
        lr.saveToVtkFile(args.vtk_file)

if __name__ == '__main__':
    main()