import numpy
import defopt
import mint
import netCDF4
import vtk

A_EARTH = 6371e3 # metres

class LFRiz(object):

    def __init__(self, infile, inmesh, pts, ndays, tindex, level):
        
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
        for icell in range(self.srcGrid.getNumberOfCells()):
            for ie in range(4):
                edgeIndex, _ = self.srcGrid.getEdgeId(icell, ie)
                n0, n1 = self.srcGrid.getNodeIds(icell, ie)
                lon0, lat0 = lon[n0], lat[n0]
                lon1, lat1 = lon[n1], lat[n1]
                latmid = 0.5*(lat0 + lat1)
                dlon, dlat = lon1 - lon0, lat1 - lat0
                self.influxes[icell, ie] = A_EARTH*(u1[edgeIndex]*dlat - u2[edgeIndex]*numpy.cos(latmid*numpy.pi/180.)*dlon)

        # create structured grids with nodal flux data attached
        numRibbons, nt = pts.shape[0] - 1, ndays + 1
        assert(numRibbons >= 1)
        self.ribbons = [numpy.zeros((2, nt), numpy.float64) for i in range(numRibbons)]
        self.fields = [numpy.zeros((2, nt), numpy.float64) for i in range(numRibbons)]


    def advect(self):
        # create a vector field interpolator
        self.vectInterp = mint.VectorInterp()
        self.vectInterp.setGrid(self.srcGrid)
        self.vectInterp.buildLocator(numCellsPerBucket=100, periodX=360.)


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
    xyz[:, :1] = xy0[:, :1]
    lfr = LFRiz(infile, inmesh, xyz, ndays, tindex, level)
    lfr.advect()
    lfr.save(outfile)

if __name__ == '__main__':
    defopt.run(main)
