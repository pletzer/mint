import netCDF4
import numpy
import vtk
from reader_base import ReaderBase

def getData(nc, mname, loc):
    """
    Find all the variables that are attached to a mesh name
    @param nc netCDF4 file handle
    @param mname mesh name
    @param loc "node", "edge", or "face"
    @return {vname: var}
    """
    res = {}
    for vname, var in nc.variables.items():
        meshName = getattr(var, 'mesh', '')
        location = getattr(var, 'location', '')
        if meshName == mname and location == loc:
            res[vname] = var
    return res


class UgridReader(ReaderBase):


    def __init__(self, filename, regularization=True):
        """
        Constructor
        @param filename UGRID file 
        """

        super(UgridReader, self).__init__()
        
        # filename:meshname
        fname, mname = filename.split(':')
        nc = netCDF4.Dataset(fname, 'r')

        mesh = nc.variables[mname][:]
        lonname, latname = nc.variables[mname].node_coordinates.split()

        # read the coordinates
        lons = nc.variables[lonname]
        lats = nc.variables[latname]

        # read the face to node connectivity
        f2nname = nc.variables[mname].face_node_connectivity
        f2ename = nc.variables[mname].face_edge_connectivity
        e2nname = nc.variables[mname].edge_node_connectivity

        f2n = nc.variables[f2nname][:]
        f2n -= nc.variables[f2nname].start_index

        f2e = nc.variables[f2ename][:]
        f2e -= nc.variables[f2ename].start_index

        e2n = nc.variables[e2nname][:]
        e2n -= nc.variables[e2nname].start_index

        ncells = f2n.shape[0]

        # construct the unstructured grid as a collection of 
        # 2D cells. Each cell has its own coordinates. Make
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

            i00, i10, i11, i01 = f2n[icell, :]

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
            if regularization and numpy.any(numpy.fabs(alts < quarterPeriodicity) < 1.e-6):
                # there is a latitude at the pole. The longitude is not well 
                # defined in this case - we can set it to any value. For 
                # esthetical reason it's good to set it to the average 
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

        # attach nodal data to the grid
        for vname, var in getData(nc, mname, 'node').items():
            # read
            data = var[:]
            nComps = 1
            if len(data.shape) > 1:
                nComps = data.shape[1]
            else:
                data = data.reshape((-1, 1))

            nCompsVec = 1
            if nComps > 1:
                # it's a vector and we need 3 components
                nCompsVec = 3

            # 4 nodes per cell, may be vector
            cData = numpy.zeros((ncells, 4, nCompsVec), numpy.float64)

            for icell in range(ncells): 
                i00, i10, i11, i01 = f2n[icell, :]
                d00, d10, d11, d01 = data[i00, :], data[i10, :], data[i11, :], data[i01, :]
                cData[icell, 0, :nComps] = d00
                cData[icell, 1, :nComps] = d10
                cData[icell, 2, :nComps] = d11
                cData[icell, 3, :nComps] = d01

            if nComps == 1:
                # scalar, remove the dimension
                cData = cData.reshape((ncells, 4,))

            print('setting point field "{}"'.format(vname))
            self.setPointField(vname, cData)

        # set edge fields
        for vname, var in getData(nc, mname, 'edge').items():

            # read
            data = var[:]

            cData = numpy.zeros((ncells, 4), numpy.float64)

            # iterate over faces
            for icell in range(ncells):

                # get the point Ids
                ptIds = f2n[icell, :]

                # iterate over the edges of the face
                for edgeId in f2e[icell, :]:

                    # get the start/end point Ids
                    ptId0, ptId1 = e2n[edgeId, :]

                    # get the face Id (0 <= ie < 4) and the sign
                    ie = -1
                    sign = 0.0

                    for i in range(4):

                        i0 = i
                        i1 = (i + 1) % 4
                        if i // 2 >= 1:
                            # last two edges
                            i0 = (i + 1) % 4
                            i1 = i

                        #       2
                        #   3--->----2
                        #   |        |
                        # 3 ^        ^ 1
                        #   |        |
                        #   0--->----1
                        #       0      

                        if ptIds[i0] == ptId0 and ptIds[i1] == ptId1:
                            # edge in netcdf file has the same orientation
                            ie = i0
                            sign = 1.0
                            break

                        elif ptIds[i0] == ptId1 and ptIds[i1] == ptId0:
                            # opposite orientation
                            ie = i1
                            sign = -1.0
                            break

                    cData[icell, ie] = sign * data[edgeId]

            print('setting edge field "{}"'.format(vname))
            self.setEdgeField(vname, cData)



        nc.close()

            



###############################################################################

def main():
    import argparse
    from numpy import pi, sin, cos, exp, heaviside, power

    parser = argparse.ArgumentParser(description='Read ugrid file')
    parser.add_argument('-i', dest='input', default='', help='Specify UGRID input netCDF file in the form FILENAME:MESHNAME')
    parser.add_argument('-r', dest='regularization', action='store_false', help='Disable pole regularization (required for latlon grid)')
    parser.add_argument('-V', dest='vtk_file', default='lonlat.vtk', help='Save grid in VTK file')
    parser.add_argument('-stream', dest='streamFunc', default='', 
                        help='Stream function as a function of x (longitude) and y (latitude)')
    parser.add_argument('-b', dest='binary', action='store_true', help='Write binary file')
    parser.add_argument('-u', dest='uFunc', default='', 
                        help='u vector component function of x (longitude) and y (latitude) on points')
    parser.add_argument('-v', dest='vFunc', default='', 
                        help='v vector component function of x (longitude) and y (latitude) on points')
    parser.add_argument('-p', dest='print', action='store_true', 
                        help='Print point coordinates cell by cell')



   
    args = parser.parse_args()

    reader = UgridReader(filename=args.input, regularization=args.regularization)

    if args.print:
        reader.printCellPoints()

    # compute the edge velocity if user provides the stream function
    x, y = reader.getLonLat()

    if args.uFunc and args.vFunc:
        uData = eval(args.uFunc)
        vData = eval(args.vFunc)
        reader.setPointVectorField('velocity_vector', uData, vData)

    if args.streamFunc:

        streamData = eval(args.streamFunc)
        print('setting point {} data to "streamFunction"'.format(streamData.shape))
        reader.setPointField('streamFunction', streamData)

        edgeVel = reader.getEdgeFieldFromStreamData(streamData)
        print('setting edge {} data to "edgeIntegratedVelocity"'.format(edgeVel.shape))
        reader.setEdgeField('edgeIntegratedVelocity', edgeVel)

        loopIntegrals = reader.getLoopIntegralsFromStreamData(streamData)
        reader.setLoopIntegrals('cell_loop_integrals', loopIntegrals)

    if args.vtk_file:
        reader.saveToVtkFile(args.vtk_file, binary=args.binary)


if __name__ == '__main__':
    main()
