import netCDF4
import numpy
import vtk
from reader_base import ReaderBase

class UgridReader(ReaderBase):


    def __init__(self, filename):
        """
        Constructor
        @param filename UGRID file 
        """

        super(UgridReader, self).__init__()
        
        # read UGRID file
        nc = netCDF4.Dataset(filename, 'r')

        lats, lons = None, None
        connectivity = None
        for varname in nc.variables:
            var = nc.variables[varname]
            if hasattr(var, 'cf_role') and var.cf_role == 'face_node_connectivity':
                connectivity = var[:]
            elif hasattr(var, 'standard_name'):
                if var.standard_name == 'longitude' and hasattr(var, 'long_name') and var.long_name.find('node') >= 0:
                    lons = var[:]
                    #print('found longitude: {}'.format(varname))
                elif var.standard_name == 'latitude' and hasattr(var, 'long_name') and var.long_name.find('node') >= 0:
                    lats = var[:]
                    #print('found latitude: {}'.format(varname))

        ncells = connectivity.shape[0]

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
        halfPeriodicity = self.PERIODICITY_LENGTH/2.
        quarterPeriodicity = self.PERIODICITY_LENGTH/4.
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
                index10 = numpy.argmin([abs(lon10 - self.PERIODICITY_LENGTH - lon00), abs(lon10 - lon00), abs(lon10 + self.PERIODICITY_LENGTH - lon00)])
                index11 = numpy.argmin([abs(lon11 - self.PERIODICITY_LENGTH - lon00), abs(lon11 - lon00), abs(lon11 + self.PERIODICITY_LENGTH - lon00)])
                index01 = numpy.argmin([abs(lon01 - self.PERIODICITY_LENGTH - lon00), abs(lon01 - lon00), abs(lon01 + self.PERIODICITY_LENGTH - lon00)])

                lon10 += (index10 - 1) * self.PERIODICITY_LENGTH
                lon11 += (index11 - 1) * self.PERIODICITY_LENGTH
                lon01 += (index01 - 1) * self.PERIODICITY_LENGTH

            lts = numpy.array([lat00, lat10, lat11, lat01])
            lns = numpy.array([lon00, lon10, lon11, lon01])
            alts = numpy.fabs(lts)
            if numpy.any(alts[:] == quarterPeriodicity):
                # there is a latitude at the pole. The longitude is not well 
                # defined in this case - we can set it to any value. For 
                # esthetical reason it't good to set it to the average 
                # of the longitudes
                i = numpy.argmax(alts - quarterPeriodicity)
                # compute the average lon value, excluding this one
                # and set lns[index] to that value
                avgLon = numpy.sum([lns[(i + 1) % 4], lns[(i + 2) % 4], lns[(i + 3) % 4]]) / 3.
                lns[i] = avgLon
                lon00, lon10, lon11, lon01 = lns

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


###############################################################################

def main():
    import argparse
    from numpy import pi, sin, cos, exp

    parser = argparse.ArgumentParser(description='Read ugrid file')
    parser.add_argument('-i', dest='input', default='', help='Specify UGRID input netCDF file')
    parser.add_argument('-p', dest='padding', type=int, default=0, 
                              help='Specify by how much the grid should be padded on the high lon side')
    parser.add_argument('-V', dest='vtk_file', default='lonlat.vtk', help='Save grid in VTK file')
    parser.add_argument('-stream', dest='streamFunc', default='', 
                        help='Stream function as a function of x (longitude) and y (latitude)')
    parser.add_argument('-b', dest='binary', action='store_true', help='Write binary file')
    parser.add_argument('-u', dest='uFunc', default='', 
                        help='u vector component function of x (longitude) and y (latitude)')
    parser.add_argument('-v', dest='vFunc', default='', 
                        help='v vector component function of x (longitude) and y (latitude)')


   
    args = parser.parse_args()

    reader = UgridReader(filename=args.input)

    # compute the edge velocity if user provides the stream function
    x, y = reader.getLonLat()

    if args.uFunc and args.vFunc:
        uData = eval(args.uFunc)
        vData = eval(args.vFunc)
        reader.setPointVectorField('velocity_vector', uData, vData)

    if args.streamFunc:

        streamData = eval(args.streamFunc)
        reader.setPointField('stream_function', streamData)

        edgeVel = reader.getEdgeFieldFromStreamData(streamData)
        reader.setEdgeField('edge_integrated_velocity', edgeVel)

        loopIntegrals = reader.getLoopIntegralsFromStreamData(streamData)
        reader.setLoopIntegrals('cell_loop_integrals', loopIntegrals)

    if args.vtk_file:
        reader.saveToVtkFile(args.vtk_file, binary=args.binary)


if __name__ == '__main__':
    main()
