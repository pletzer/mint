import numpy
import netCDF4
import argparse

class LatLon:

    def __init__(self):
        """
        Constructor
        no arguments
        """

        # earth radius
        self.a = 6371.e3

        self.fillValue = -1.073742e+09

    def setEarthRadius(self, a):
        """
        Set the earth's radius in metre
        @param a radius
        """
        self.a = a

    def setNumberOfLonCells(self, numLons):
        """
        Set the number of longitudinal cells
        @param numLons number of cells
        """
        self.numLons0 = numLons

    def setNumberOfLatCells(self, numLats):
        """
        Set the number of latitudinal cells
        @param numLats number of cells
        """
        self.numLats0 = numLats

    def build(self):
        """
        Build the object. Call this after setNumberOfLatCells and setNumberOfLonCells
        """

        # deltas
        self.dLat = 180.0/float(self.numLats0)
        self.dLon = 360.0/float(self.numLons0)

        # latitude/longitude array
        self.lats = numpy.linspace(-90.0, 90.0, self.numLats0 + 1)
        self.lons = numpy.linspace(0., 360.0, self.numLons0 + 1)

        # vector field components
        self.u = self.fillValue * numpy.ones((self.numLats0, self.numLons0 + 1), numpy.float64)
        self.v = self.fillValue * numpy.ones((self.numLats0 + 1, self.numLons0), numpy.float64)


    def setStreamFunction(self, expr):
        """
    	Set vector field data from stream function
    	@param expr expression of x, y for u component (x is longitude, y is latitude)
    	"""
        y, x = numpy.meshgrid(self.lats, self.lons, indexing='ij')
        psi = eval(expr)

        # divide by the edge length to get the intensive variables
        dy = self.a * self.dLat * numpy.pi/180.
        dx = self.a * self.dLon * numpy.pi/180.
        self.u[...] = (psi[1:, :] - psi[:-1, :]) / dy
        self.v[...] = (psi[:, :-1] - psi[:, 1:]) / dx


    def load(self, filename):
        """
        Load data from NetCDF file
        @param filename file name
        """
        nc = netCDF4.Dataset(filename, 'r')

        self.a = float(nc.earth_radius)

        self.numLats0 = nc.dimensions['latitude_0'][:]
        self.numLons0 = nc.dimensions['longitude_0'][:]

        self.lats = numpy.zeros((self.numLats0 + 1,), numpy.float64)
        self.lons = numpy.zeros((self.numLons0 + 1,), numpy.float64)

        self.lats[:] = nc.variables['latitude'][:]
        # we're not saving last value because of periodicity
        self.lons[:-1] = nc.variables['longitude'][:]

        self.u = nc.variables['u'][:]
        self.v = nc.variables['v'][:]

        nc.close()


    def dump(self, filename):

        nc = netCDF4.Dataset(filename, 'w')

        nc.earth_radius = str(self.a) + " metres"
        
        # UM does not store the last longitude
        longitude_dim = nc.createDimension('longitude', self.numLons0)
        latitude_dim = nc.createDimension('latitude', self.numLats0 + 1)
        longitude_0_dim = nc.createDimension('longitude_0', self.numLons0)
        latitude_0_dim = nc.createDimension('latitude_0', self.numLats0)

        longitude_var = nc.createVariable('longitude', 'f8', ('longitude',))
        latitude_var = nc.createVariable('latitude', 'f8', ('latitude',))
        longitude_0_var = nc.createVariable('longitude_0', 'f8', ('longitude_0',))
        latitude_0_var = nc.createVariable('latitude_0', 'f8', ('latitude_0',))

        latitude_var.axis = "Y"
        latitude_var.units = "degrees_north"
        latitude_var.standard_name = "latitude"
        longitude_var.axis = "X"
        longitude_var.units = "degrees_east"
        longitude_var.standard_name = "longitude"

        latitude_0_var.axis = "Y"
        latitude_0_var.units = "degrees_north"
        latitude_0_var.standard_name = "latitude"
        longitude_0_var.axis = "X"
        longitude_0_var.units = "degrees_east"
        longitude_0_var.standard_name = "longitude"        

        # write the data
        latitude_var[:] = self.lats
        latitude_0_var[:] = self.lats[:-1] + 0.5*self.dLat

        longitude_var[:] = self.lons[:-1] # skip the last value
        longitude_0_var[:] = self.lons[:-1] + 0.5*self.dLon

        u_var = nc.createVariable('u', 'f8', ('latitude_0', 'longitude',), fill_value=self.fillValue)
        u_var.standard_name = 'eastward_wind'
        u_var.units = 'm s-1'
        u_var.grid_mapping = 'latitude_longitude'
        u_var[:] = self.u[:, :-1]

        v_var = nc.createVariable('v', 'f8', ('latitude', 'longitude_0',), fill_value=self.fillValue)
        v_var.standard_name = 'northward_wind'
        v_var.units = 'm s-1'
        v_var.grid_mapping = 'latitude_longitude'
        v_var[:] = self.v[:, :]

        nc.close()


#############################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-nlat', dest='nlat', type=int, default=2, help='Number of latitude cells')
    parser.add_argument('-nlon', dest='nlon', type=int, default=4, help='Number of longitude cells')
    parser.add_argument('-strmfc', dest='strmfc', type=str, default='x/360.', help='Stream function of x (longitude)  and y (latitude)')
    parser.add_argument('-o', dest='output', type=str, default='lalon.nc', help='Output file')

    args = parser.parse_args()
    numLat, numLon = args.nlat, args.nlon
    streamFunc = args.strmfc
    outputFile = args.output

    ll = LatLon()
    ll.setNumberOfLatCells(numLat)
    ll.setNumberOfLonCells(numLon)
    ll.build()
    ll.setStreamFunction(streamFunc)
    ll.dump(outputFile)

if __name__ == '__main__':
    main()

