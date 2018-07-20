import netCDF4
import numpy
import vtk
from reader_base import ReaderBase

class LatLonReader(ReaderBase):


    def __init__(self, filename, padding=0):
        """
        Constructor
        @param filename UM netCDF file
        @param padding number of extra cells to add on the high end of longitudes
        @note padding add extra cells on the high end of longitudes
        """

        super(LatLonReader, self).__init__()
        
        # read UGRID file
        nc = netCDF4.Dataset(filename, 'r')

        lon_units = ''
        lat_units = ''

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
                    lons_units = var.units
                elif var.standard_name == 'latitude':
                    if varname.find('_0') >= 0:
                        lats_0 = var[:]
                    else:
                        lats = var[:]
                    lats_units = var.units

        ncells_lat, ncells_lon = len(lats_0), len(lons_0)
        ncells = ncells_lat * (ncells_lon + padding)

        # construct the unstructured grid as a collection of 
        # 2D cells

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

        periodicity_length = 360. # in deg
        icell = 0
        for j0 in range(ncells_lat):

            j1 = j0 + 1

            for i in range(ncells_lon + padding):

                i0 = (i + 0) % ncells_lon
                i1 = (i + 1) % ncells_lon
                offset0 = periodicity_length * ((i + 0) // ncells_lon)
                offset1 = periodicity_length * ((i + 1) // ncells_lon)

                lon00, lat00 = lons[i0] + offset0, lats[j0]
                lon10, lat10 = lons[i1] + offset1, lats[j0]
                lon11, lat11 = lons[i1] + offset1, lats[j1]
                lon01, lat01 = lons[i0] + offset0, lats[j1]

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


###############################################################################

def main():
    import argparse
    from numpy import pi, cos, sin, exp

    parser = argparse.ArgumentParser(description='Read ugrid file')
    parser.add_argument('-i', dest='input', default='ll.nc', help='Specify UM input netCDF file')
    parser.add_argument('-p', dest='padding', type=int, default=0, 
                              help='Specify by how much the grid should be padded on the high lon side')
    parser.add_argument('-V', dest='vtk_file', default='lonlat.vtk', help='Save grid in VTK file')
    parser.add_argument('-b', dest='binary', action='store_true', help='Write binary file')
    parser.add_argument('-stream', dest='streamFunc', default='x', 
                        help='Stream function as a function of x (longitude in rad) and y (latitude in rad)')

   
    args = parser.parse_args()

    reader = LatLonReader(filename=args.input, padding=args.padding)

    if args.streamFunc:
        # compute the edge velocity if user provides the stream function
        x, y = reader.getLonLat()
        streamData = eval(args.streamFunc)
        edgeVel = reader.getEdgeFieldFromStreamData(streamData)
        reader.setEdgeField('edge_integrated_velocity', edgeVel)
        loopIntegrals = reader.getLoopIntegralsFromStreamData(streamData)
        reader.setLoopIntegrals('cell_loop_integrals', loopIntegrals)

    if args.vtk_file:
        reader.saveToVtkFile(args.vtk_file, binary=args.binary)

if __name__ == '__main__':
    main()
