import numpy
import defopt
import mint
import netCDF4
import vtk
from scipy.integrate import odeint

A_EARTH = 6371.e3 # metres
DEG2RAD = numpy.pi/180.

class LFRiz(object):

    def __init__(self, infile, inmesh, pts, ndays, tindex, level):

        self.ndays = ndays
        self.nt = ndays + 1
        print(f'num times = {self.nt}')
        self.pts0 = pts

        self.numPoints = self.pts0.shape[0]

        # read the data
        self.srcGrid = mint.Grid()
        # cubed-sphere
        self.srcGrid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1)
        self.srcGrid.loadFromUgrid2D(infile + '$' + inmesh)

        # read the u, v components for each unique edge, in m/s
        nc = netCDF4.Dataset(infile)
        u1 = nc.variables['u1'][tindex, level, :]
        u2 = nc.variables['u2'][tindex, level, :]

        # read the coordinates
        xname, yname = nc.variables[inmesh].node_coordinates.split(' ')
        lon = nc.variables[xname][:]
        lat = nc.variables[yname][:]

        self.influxes = numpy.zeros((self.srcGrid.getNumberOfCells(), 4), numpy.float64)

        coeff = 3600 * 24 / A_EARTH

        # iterate
        for icell in range(self.srcGrid.getNumberOfCells()):
            for ie in range(4):

                edgeIndex, edgeSign = self.srcGrid.getEdgeId(icell, ie)
                n0, n1 = self.srcGrid.getNodeIds(icell, ie)

                lon0, lat0 = lon[n0], lat[n0]
                lon1, lat1 = lon[n1], lat[n1]

                latmid_rad = 0.5*(lat0 + lat1) * DEG2RAD
                dlon, dlat = lon1 - lon0, lat1 - lat0

                # fluxes in deg^2/day (approximation)
                self.influxes[icell, ie] = coeff *  (u1[edgeIndex] * dlat / numpy.cos(latmid_rad) - u2[edgeIndex] * dlon)
                
                print(f'cell {icell} edge {ie} edgeSign={edgeSign} u1,u2 = {u1[edgeIndex]},{u2[edgeIndex]} (m/s) p0={lon0},{lat0} p1={lon1},{lat1} dlon,dlat={dlon},{dlat} (deg) flux = {self.influxes[icell, ie]} (deg^2/day)')


        # compute the fluxes in deg^2/day across each edge        # from m/s to rad/day
        coeff = 3600 * 24 / A_EARTH # sec -> day

        # iterate
        for icell in range(self.srcGrid.getNumberOfCells()):
            for ie in range(4):

                edgeIndex, _ = self.srcGrid.getEdgeId(icell, ie)
                n0, n1 = self.srcGrid.getNodeIds(icell, ie)

                lon0, lat0 = lon[n0], lat[n0]
                lon1, lat1 = lon[n1], lat[n1]
                latmid = 0.5*(lat0 + lat1)
                dlon, dlat = lon1 - lon0, lat1 - lat0

                # convert to rads
                dlon *= DEG2RAD
                dlat *= DEG2RAD
                latmid *= DEG2RAD
                # fluxes in rad^2/day (approximation)
                self.influxes[icell, ie] = coeff*( u1[edgeIndex] * dlat / numpy.cos(latmid) - u2[edgeIndex] * dlon)
                #print(f'cell {icell} edge {ie} u1, u2 = {u1[edgeIndex]}, {u2[edgeIndex]} dlon, dlat = {dlon}, {dlat} [rad] flux = {self.influxes[icell, ie]} ')


    def advect(self):

        # create a vector field interpolator
        vi = mint.VectorInterp()
        vi.setGrid(self.srcGrid)
        vi.buildLocator(numCellsPerBucket=100, periodX=360.)

        def tendency(xyz, t):
            p = xyz.reshape((self.numPoints, 3))
            vi.findPoints(p)
            vect = vi.getFaceVectors(self.influxes, placement=0)
            print(f'vect = {vect}')
            return vect.reshape((self.numPoints*3,))

        # time step
        dt = 1. # one day

        # all the time values for which we seek the advected positions
        tvals = [i*dt for i in range(self.nt)]
        print(f'time values: {tvals}')

        # intial positions as a flat vector
        xyz0 = self.pts0.reshape((self.numPoints*3,))
        print(f'initial positions = {xyz0}')

        # advance the positions
        sol = odeint(tendency, xyz0, tvals, rtol=1.e-12, atol=1.e-12)

        # save the positions, column major
        self.points = sol.reshape((self.nt, self.numPoints, 3))
        print(f'advected points:')
        print(f'{self.points}')


    def build(self):

        pli = mint.PolylineIntegral() # to compute fluxes

        numRibbons = self.numPoints - 1
        assert(numRibbons >= 1) # need at least two points

        self.vfluxData = vtk.vtkDoubleArray()
        self.vfluxData.SetName('flux')
        self.vfluxData.SetNumberOfComponents(1)
        self.vfluxData.SetNumberOfTuples(numRibbons * 2 * self.nt)

        self.vpointData = vtk.vtkDoubleArray()
        self.vpointData.SetNumberOfComponents(3)
        self.vpointData.SetNumberOfTuples(numRibbons * 2 * self.nt)

        self.vpoints = vtk.vtkPoints()
        self.vpoints.SetData(self.vpointData)

        self.vgrid = vtk.vtkUnstructuredGrid()
        self.vgrid.SetPoints(self.vpoints)
        self.vgrid.Allocate(1, 1)

        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(4)

        for i in range(numRibbons):

            # set the nodal values
            for j in range(self.nt):

                index0 = 2*(self.nt*i + j)
                index1 = index0 + 1

                xyz0 = self.points[j, i + 0, :]
                xyz1 = self.points[j, i + 1, :]
                self.vpointData.SetTuple(index0, xyz0)
                self.vpointData.SetTuple(index1, xyz1)

                tpts = numpy.ascontiguousarray( self.points[j, i:i+2, :] )
                pli.build(self.srcGrid, tpts, counterclock=False, periodX=360.)
                flx = (pli.getIntegral(self.influxes),)

                # same flux values for 2 adjacent nodes
                self.vfluxData.SetTuple(index0, flx)
                self.vfluxData.SetTuple(index1, flx)

            # build the cells
            for j in range(self.ndays): # self.nt - 1

                index0 = 2*(self.nt*i + j)
                index1 = index0 + 1
                index2 = index1 + 2
                index3 = index2 - 1

                ptIds.SetId(0, index0)
                ptIds.SetId(1, index1)
                ptIds.SetId(2, index2)
                ptIds.SetId(3, index3)

                self.vgrid.InsertNextCell(vtk.VTK_QUAD, ptIds)

        self.vgrid.GetPointData().AddArray(self.vfluxData)
        print(self.vgrid)

 
    def show(self):
        pass


    def save(self, filename):

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.vgrid)
        writer.Update()

        # numRibbons = self.numPoints - 1
        # assert(numRibbons >= 1) # need at least two points

        # for i in range(numRibbons):

        #     vfluxData = vtk.vtkDoubleArray()
        #     vpointData = vtk.vtkDoubleArray()
        #     vpoints = vtk.vtkPoints()
        #     vsgrid = vtk.vtkStructuredGrid()
        #     vwriter = vtk.vtkStructuredGridWriter()

        #     # set the dimensions
        #     vfluxData.SetNumberOfTuples(self.nt * 2)
        #     vfluxData.SetNumberOfComponents(1) # scalar
        #     vpointData.SetNumberOfTuples(self.nt * 2)
        #     vpointData.SetNumberOfComponents(3)
        #     vpoints.SetNumberOfPoints(self.nt * 2)
        #     vsgrid.SetDimensions(2, self.nt, 1)
        #     fname = filename + f'_{i}.vtk'
        #     vwriter.SetFileName(fname)

        #     # vertices
        #     pts = numpy.ascontiguousarray( self.points[:, i:i+2, :] )

        #     # fluxes
        #     fluxes = numpy.zeros( (self.nt, 2), dtype=numpy.float64)


        #     pli = mint.PolylineIntegral()
        #     for j in range(self.nt):
        #         tpts = numpy.ascontiguousarray( self.points[j, i:i+2, :] )
        #         pli.build(self.srcGrid, tpts, counterclock=False, periodX=360.)
        #         fluxes[j, 0:2] = pli.getIntegral(self.influxes)


        #     # connect
        #     vfluxData.SetVoidArray(fluxes, self.nt * 2, 1)
        #     vfluxData.SetName('flux')
        #     vpointData.SetVoidArray(pts, self.nt * 2 * 3, 1)
        #     vpoints.SetData(vpointData)
        #     vsgrid.SetPoints(vpoints)
        #     vsgrid.GetPointData().AddArray(vfluxData)
        #     vwriter.SetInputData(vsgrid)
        #     vwriter.Update()

        #     print(f'wrote {fname}')






###############################################################################

def main(*,
         infile: str = '../data/lfric_diag_wind.nc', inmesh: str = 'Mesh2d',
         pts: str="(20., 40.), (50., 40.)", 
         ndays : int=10, 
         tindex : int = 0, level: int = 0, 
         outfile: str = 'lfriz.vtk'):

    """
    Generate advection grid

    param: infile netcdf file containing u1, u2 field
    param: inmesh mesh name in the netcdf file
    param: pts set of initial points
    param: ndays number of days to integrate forward
    param: tindex time index
    param: level level
    param: outfile output VTK file
    """

    xy0 = numpy.array(eval(pts))
    xyz = numpy.zeros((xy0.shape[0], 3))
    xyz[:, 0:2] = xy0[:, 0:2]
    lfr = LFRiz(infile, inmesh, xyz, ndays, tindex, level)
    lfr.advect()
    lfr.build()
    lfr.save(outfile)

if __name__ == '__main__':
    defopt.run(main)
