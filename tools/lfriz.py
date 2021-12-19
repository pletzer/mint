import numpy
import defopt
import mint
import netCDF4

class LFRiz(object):

    def __init__(self, infile, inmesh, pts, ndays, tindex, level):
        
        # read the data
        self.srcGrid = mint.Grid()
        # cubed-sphere
        self.srcGrid.setFlags(fixLonAcrossDateline=1, averageLonAtPole=1)
        self.srcGrid.loadFromUgrid2D(infile + inmesh)

        # associate fluxes to each edge
        nc = netCDF4.Dataset(infile)
        u1 = nc.variables['u1'][:]
        u2 = nc.variables['u2'][:]

        # create a vector field interpolator
        self.vectInterp = mint.VectorInterp()
        self.vectInterp.setGrid(self.srcGrid)
        self.vectInterp.buildLocator(numCellsPerBucket=100, periodX=360.)

        # create structured grids


    def advect(self):
        pass


    def show(self):
        pass


    def save(self, filename):
        pass

###############################################################################

def main(*,
         infile: str = '../data/lfric_diag_wind.nc', inmesh: str = '$Mesh2d',
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
