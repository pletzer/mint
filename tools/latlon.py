import numpy
import netCDF4
import argparse

"""
A class to generate a uniform lat-lon grid conforming to the output of UM
"""

class LatLon:

    def __init__(self):
        """
        Constructor
        no arguments
        """
        self.fillValue = -1.073742e+09


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


    def load(self, filename):
        """
        Load data from NetCDF file
        @param filename file name
        """
        nc = netCDF4.Dataset(filename, 'r')

        self.numLats0 = nc.dimensions['latitude_0'].size
        self.numLons0 = nc.dimensions['longitude_0'].size

        self.lats = numpy.zeros((self.numLats0 + 1,), numpy.float64)
        self.lons = numpy.zeros((self.numLons0 + 1,), numpy.float64)

        self.lats[:] = nc.variables['latitude'][:]
        # we're not saving the last value because of periodicity
        self.lons[:-1] = nc.variables['longitude'][:]

        nc.close()


    def dump(self, filename):

        nc = netCDF4.Dataset(filename, 'w')
        
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

        nc.close()


#############################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-nlat', dest='nlat', type=int, default=2, help='Number of latitude cells')
    parser.add_argument('-nlon', dest='nlon', type=int, default=4, help='Number of longitude cells')
    parser.add_argument('-o', dest='output', type=str, default='um.nc', help='Output file')

    args = parser.parse_args()
    numLat, numLon = args.nlat, args.nlon
    outputFile = args.output

    ll = LatLon()
    ll.setNumberOfLatCells(numLat)
    ll.setNumberOfLonCells(numLon)
    ll.build()
    ll.dump(outputFile)

if __name__ == '__main__':
    main()

