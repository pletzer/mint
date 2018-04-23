import numpy
import netCDF4

class LatLon:

    def __init__(self, numLats, numLons):

        # create grid
        self.numLats0 = numLats
        self.numLons0 = numLons
        self.numLats1 = numLats + 1
        self.numLons1 = numLons

        self.dLat = 180.0/float(numLats)
        self.dLon = 360.0/float(numLons)


        # latitudes
        self.lats = numpy.linspace(-90.0, 90.0, self.numLats1)

        # longitudes: last cell implicitly wraps around 
        self.lons = numpy.linspace(0., 360.0 - self.dLon, self.numLons1)


    def saveToNetcdf(self, filename, u=None, v=None):

        nc = netCDF4.Dataset(filename, 'w')
        
        longitude_dim = nc.createDimension('longitude', self.numLons1)
        latitude_dim = nc.createDimension('latitude', self.numLats1)
        longitude_0_dim = nc.createDimension('longitude_0', self.numLons0)
        latitude_0_dim = nc.createDimension('latitude_0', self.numLats0)

        longitude_var = nc.createVariable('longitude', 'f4', ('longitude',))
        latitude_var = nc.createVariable('latitude', 'f4', ('latitude',))
        longitude_0_var = nc.createVariable('longitude_0', 'f4', ('longitude_0',))
        latitude_0_var = nc.createVariable('latitude_0', 'f4', ('latitude_0',))

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
        longitude_var[:] = self.lons
        latitude_0_var[:] = self.lats[:-1] + 0.5*self.dLat
        longitude_0_var[:] = self.lons + 0.5*self.dLon

        if u:
        	u_var = nc.createVariable('u', 'f4', ('latitude_0', 'longitude',))
        	u_var._FillValue = -1.073742e+09
        	u_var.standard_name = 'eastward_wind'
        	u_var.units = 'm s-1'
        	u_var.grid_mapping = 'latitude_longitude'
        	u_var[:] = u

        if v:
        	v_var = nc.createVariable('v', 'f4', ('latitude', 'longitude_0',))
        	v_var._FillValue = -1.073742e+09
        	v_var.standard_name = 'northward_wind'
        	v_var.units = 'm s-1'
        	v_var.grid_mapping = 'latitude_longitude'
        	v_var[:] = v

    	nc.close()


#############################################################################
def test():
    numLat, numLon = 2, 4
    ll = LatLon(numLat, numLon)
    ll.saveToNetcdf('ll.nc')

if __name__ == '__main__':
    test()

