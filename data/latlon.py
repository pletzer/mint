import numpy
import netCDF4

class LatLon:

    def __init__(self, numLats, numLons):

        # create grid
        self.numLats0 = numLats
        self.numLons0 = numLons
        self.numLats1 = numLats + 1
        self.numLons1 = numLons + 1

        self.dLat = 180.0/float(numLats)
        self.dLon = 360.0/float(numLons)


        # latitudes
        self.lats = numpy.linspace(-90.0, 90.0, self.numLats1)

        # longitudes: last cell implicitly wraps around 
        self.lons = numpy.linspace(0., 360.0, self.numLons1)

        self.fillValue = -1.073742e+09

        self.uDataExtensive = self.fillValue * numpy.ones((self.numLats0, self.numLons1), numpy.float64)
        self.vDataExtensive = self.fillValue * numpy.ones((self.numLats1, self.numLons0), numpy.float64)

        self.uDataIntensive = self.fillValue * numpy.ones((self.numLats0, self.numLons1), numpy.float64)
        self.vDataIntensive = self.fillValue * numpy.ones((self.numLats1, self.numLons0), numpy.float64)

        # earth radius
        self.a = 6371.e3

    def setEarthRadius(self, a):
        self.a = a


    def setStreamFunction(self, expr):
        """
    	Set vector field data from stream function
    	@param expr expression of x, y for u component (x is longitude, y is latitude)
    	"""
        y, x = numpy.meshgrid(self.lats, self.lons, indexing='ij')
        psi = eval(expr)

        # line integrated values
        self.uDataExtensive[...] = psi[1:, :] - psi[:-1, :]
        self.vDataExtensive[...] = psi[:, :-1] - psi[:, 1:]

        # divide by the edge length to get the intensive variables
        dy = self.a * self.dLat * numpy.pi/180.
        dx = self.a * self.dLon * numpy.pi/180.
        self.uDataIntensive[...] = self.uDataExtensive[...] / dy
        self.vDataIntensive[...] = self.vDataExtensive[...] / dx


    def saveToNetcdf(self, filename):

        nc = netCDF4.Dataset(filename, 'w')
        
        # UM does not store the last longitude
        longitude_dim = nc.createDimension('longitude', self.numLons1 - 1)
        latitude_dim = nc.createDimension('latitude', self.numLats1)
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
        longitude_var[:] = self.lons[:-1] # skip the last value
        latitude_0_var[:] = self.lats[:-1] + 0.5*self.dLat
        longitude_0_var[:] = self.lons[:-1] + 0.5*self.dLon

        u_var = nc.createVariable('u', 'f8', ('latitude_0', 'longitude',), fill_value=self.fillValue)
        u_var.standard_name = 'eastward_wind'
        u_var.units = 'm s-1'
        u_var.grid_mapping = 'latitude_longitude'
        u_var[:] = self.uDataIntensive[:, :-1]

        v_var = nc.createVariable('v', 'f8', ('latitude', 'longitude_0',), fill_value=self.fillValue)
        v_var.standard_name = 'northward_wind'
        v_var.units = 'm s-1'
        v_var.grid_mapping = 'latitude_longitude'
        v_var[:] = self.vDataIntensive[:, :]

    	nc.close()


#############################################################################
def test():
    numLat, numLon = 2, 4
    ll = LatLon(numLat, numLon)
    ll.setStreamFunction('x')
    ll.saveToNetcdf('latlon.nc')

if __name__ == '__main__':
    test()

