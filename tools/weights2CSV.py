import pandas
import netCDF4
import argparse

parser = argparse.ArgumentParser(description='Convert weights from netCDF to CSV')
parser.add_argument('-i', dest='weights_nc_file', default='',
                    help='Specify input NetCDF file')
parser.add_argument('-o', dest='weights_cvs_file', default='',
                    help='Specify output CSV file')
args = parser.parse_args()


nc = netCDF4.Dataset(args.weights_nc_file)
data = {}
for vname in nc.variables:
    # skip anything that has param in its name
    if vname.find('param') < 0:
        # not a parameter
        data[vname] = nc.variables[vname][:].data

csv_data = pandas.DataFrame(data)
if args.weights_cvs_file:
    csv_data.to_csv(args.weights_cvs_file)
else:
    print(csv_data)



