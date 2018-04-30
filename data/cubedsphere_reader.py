import netCDF4
import numpy
import vtk

class CubedsphereReader:

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

        points = vtk.vtkPoints()
        points.SetNumberOfPoints(4 * ncells)

        grid = vtk.vtkUnstructuredGrid()
        grid.Allocate(1, 1)
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

            k = 4*icell

            # storing coords as lon, lat, 0
            points.InsertPoint(k + 0, lon00, lat00, 0.)
            points.InsertPoint(k + 1, lon10, lat10, 0.)
            points.InsertPoint(k + 2, lon11, lat11, 0.)
            points.InsertPoint(k + 3, lon01, lat01, 0.)

            ptIds.SetId(0, k + 0)
            ptIds.SetId(1, k + 1)
            ptIds.SetId(2, k + 2)
            ptIds.SetId(3, k + 3)
            grid.InsertNextCell(vtk.VTK_QUAD, ptIds)

        grid.SetPoints(points)

        # add a cell locator
        loc = vtk.vtkCellLocator()
        loc.SetDataSet(grid)
        loc.BuildLocator()

        # store
        self.vtk = {
            'points': points,
            'grid': grid,
            'locator': loc,
        }


        # for finding intersections
        self.p0 = numpy.zeros((3,), numpy.float64)
        self.p1 = numpy.zeros((3,), numpy.float64)
        self.point = numpy.zeros((3,), numpy.float64)
        self.pcoords = numpy.zeros((3,), numpy.float64)
        self.t = vtk.mutable(-1.)
        self.subId = vtk.mutable(-1)
        self.cellId = vtk.mutable(-1)
        self.weights = numpy.zeros((4,), numpy.float64)
        self.cell = vtk.vtkGenericCell()

        
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


    def findCell(self, lonlat, tol=1.e-10):
        """
        Find cell containing point 
        @param lonlat target point
        @param tol tolerance
        @return cell Id
        """
        self.point[:2] = lonlat
        cId = self.vtk['locator'].FindCell(self.point, tol, self.cell, self.pcoords, self.weights)
        return cId


    def computeLineIntersectionPoints(self, lonlat0, lonlat1, tol=1.e-10):
        """
        Compute all the intersection points with given line
        @param lonlat0 starting point of the line
        @param lonlat1 end point of the line
        @return list [(cellId0, xi0, eta0, t0), (cellId1, xi1, eta1, t1), ...],
                     where t0, t1 are the parametric coordinates along the line
        """
        res = []
        tStart = 0.0
        loc = self.vtk['locator']
        cell = vtk.vtkGenericCell()

        # need to invert the order!!!!!
        self.p0[:2] = lonlat0[1], lonlat0[0]
        self.p1[:2] = lonlat1[1], lonlat1[0]

        # always add the starting point
        cId = loc.FindCell(self.p0, tol, cell, self.pcoords, self.weights)
        if cId >= 0:
            res.append( (cId, self.pcoords[0], self.pcoords[1], 0.0)  )
        else:
            print('Warning: starting point {} not found!'.format(lonlat0))

        # find all intersection points in between
        found = True
        while found:
            found = loc.IntersectWithLine(self.p0, self.p1, tol, self.t, 
                                          self.point, self.pcoords, self.subId, self.cellId, cell)
            if found:

                cId = self.cellId.get()

                # correct the parametric coordinate to account for moving the starting point
                # (self.t.get() is the param coord from self.p0 to lonlat1 with self.p0
                # moving forward)
                t = tStart + (self.t.get() - tStart)/(1.0 - tStart)

                # add the contribution
                res.append( (cId, self.pcoords[0], self.pcoords[1], t) ) 

                # reset the start position
                tStart = t
                self.p0[:] = self.point + eps*(self.p1 - self.p0)

        # always add the endpoint
        cId = loc.FindCell(self.p1, tol, cell, self.pcoords, self.weights)
        if cId >= 0:
            res.append( (cId, self.pcoords[0], self.pcoords[1], 1.0) )
        else:
            print('Warning: end point {} not found!'.format(lonlat1))

        return res





###############################################################################

def main():
    import argparse
    from math import pi

    parser = argparse.ArgumentParser(description='Read cubedsphere file')
    parser.add_argument('-i', dest='input', default='', help='Specify input file')
    parser.add_argument('-V', dest='vtk_file', default='', help='Save grid in VTK file')
    parser.add_argument('-p0', dest='p0', default='0., 0.', help='Starting position for target line')
    parser.add_argument('-p1', dest='p1', default='0., 2*pi', help='End position for target line')
    parser.add_argument('-p', dest='point', default='5.5, 0.3', help='Target point (lon, lat)')    
   
    args = parser.parse_args()

    csr = CubedsphereReader(filename=args.input)

    point = eval(args.point)
    cId = csr.findCell(point)
    print 'point {} is in cell {}'.format(point, cId)


    latlon0 = eval(args.p0)
    latlon1 = eval(args.p1)
    interPoints = csr.computeLineIntersectionPoints(latlon0, latlon1)
    print interPoints

    if args.vtk_file:
        csr.saveToVtkFile(args.vtk_file)

if __name__ == '__main__':
    main()