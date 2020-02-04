# loop through files to test build_concave_hulls
import os
import sys
import netCDF4
import logging
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils


logger = logging.getLogger()

fh = logging.FileHandler('build_concave_hulls_from_grids.log')
fh.setLevel(logging.DEBUG)
logger.addHandler(fh)

input_dir = sys.argv[1]
num_of_files_processed = 0
num_of_files_failed = 0

for filename in os.listdir(input_dir):
    try:
        ds = netCDF4.Dataset(filename, 'r')
        ngu = NetCDFGridUtils(ds)
        shapely_shape = ngu.get_concave_hull()
        logger.debug(shapely_shape.area)
        logger.debug(shapely_shape.wkt)
        num_of_files_processed = num_of_files_processed + 1
    except Exception as e:
        logger.debug("error on file: {}".format(filename))
        logger.debug(e)
        num_of_files_failed = num_of_files_failed + 1


logger.debug("Number of files proccessd: {}".format(num_of_files_processed))
logger.debug("Number of files failed: {}".format(num_of_files_failed))
# #print(ds)
#
# #print(ds.variables)
# #values = ds['Band1']
# # latlong_tup = (ds.dimensions['lat'], ds.dimensions['lon'])
# # ds.createVariable('fake_values', 'f', latlong_tup)
# # print("DIMENSIONS")
# # print(ds.dimensions)
# # print("LAt size")
# # print(len(ds.dimensions['lat']))
# # lat_size = ds.dimensions['lat'].size
# # lon_size = ds.dimensions['lon'].size
# # print(lat_size)
# # make list of the correction dimensions and sizes
# # ds['fake_values'][:] = numpy.zeros(shape=(lat_size, lon_size))
#
#
#
# lon = ds['lon']
# lat = ds['lat']






