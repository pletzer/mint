import numpy
import defopt
import mint
import netCDF4
import vtk
from scipy.integrate import odeint

A_EARTH = 6371e3 # metres


class LFRiz(object):

    def __init__(self, infile, inmesh, pts, ndays, tindex, level):

        self.ndays = ndays
        self.nt = ndays + 1
        self.pts0 = pts

        self.numPoints = self.pts0.shape[0]

        # read the data
        self.srcGrid = mint.Grid()
        # cubed-sphere
        self.srcGrid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1)
        self.srcGrid.loadFromUgrid2D(infile + '$' + inmesh)

        # read the u, v components, in m/s
        nc = netCDF4.Dataset(infile)
        u1 = nc.variables['u1'][tindex, level, :]
        u2 = nc.variables['u2'][tindex, level, :]

        # read the edge to node connectivity
        e2nname = nc.variables[inmesh].edge_node_connectivity
        e2n = nc.variables[e2nname][:]

        # read the coordinates
        xname, yname = nc.variables[inmesh].node_coordinates.split(' ')
        lon = nc.variables[xname][:]
        lat = nc.variables[yname][:]

        # compute the fluxes in m^2/s across each edge
        uvec = numpy.zeros((2,), numpy.float64)
        self.influxes = numpy.zeros((self.srcGrid.getNumberOfCells(), 4), numpy.float64)
        # from m/s to rad/day
        coeff = numpy.pi * 3600 * 24 / (180. * A_EARTH)
        deg2rad = numpy.pi/180.
        for icell in range(self.srcGrid.getNumberOfCells()):
            for ie in range(4):
                edgeIndex, _ = self.srcGrid.getEdgeId(icell, ie)
                n0, n1 = self.srcGrid.getNodeIds(icell, ie)
                lon0, lat0 = lon[n0], lat[n0]
                lon1, lat1 = lon[n1], lat[n1]
                latmid = 0.5*(lat0 + lat1)
                dlon, dlat = lon1 - lon0, lat1 - lat0
                # convert to rads
                dlon *= deg2rad
                dlat *= deg2rad
                latmid *= deg2rad
                # fluxes in rad^2/day
                self.influxes[icell, ie] = coeff*( u1[edgeIndex] * dlat / numpy.cos(latmid) - u2[edgeIndex] * dlon)


    def advect(self):

        # create a vector field interpolator
        vi = mint.VectorInterp()
        vi.setGrid(self.srcGrid)
        vi.buildLocator(numCellsPerBucket=100, periodX=360.)

        def tendency(xyz, t):
            p = xyz.reshape((self.numPoints, 3))
            vi.findPoints(p)
            vect = vi.getFaceVectors(self.influxes, placement=0)
            return vect.reshape((self.numPoints*3,))

        # time step
        dt = 1. # one day

        # all the time values for which we seek the advected positions
        tvals = [i*dt for i in range(self.nt)]

        # intial positions as a flat vector
        xyz0 = self.pts0.reshape((self.numPoints*3,))

        # advance the positions
        sol = odeint(tendency, xyz0, tvals, rtol=1.e-12, atol=1.e-12)

        # save the positions
        self.points = sol.reshape((self.nt, self.numPoints, 3))

        print(self.points)



    def show(self):
        pass


    def save(self, filename):
        pass

###############################################################################

def main(*,
         infile: str = '../data/lfric_diag_wind.nc', inmesh: str = 'Mesh2d',
         pts: str="(0., 10.), (50., 10.)", 
         ndays : float=10, 
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
    lfr.save(outfile)

if __name__ == '__main__':
    defopt.run(main)
