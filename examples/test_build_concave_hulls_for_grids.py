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


def build_concave_on_single_file(input_file):

   # file = "C:/Users/u62231/Desktop/Projects/gadds/grid_exmples/P633demg.nc"
    ds = netCDF4.Dataset(input_file)
    logger.debug(input_file)
   # value = ds['Band1']._FillValue
    ngu = NetCDFGridUtils(ds)
    shapely_shape = ngu.get_concave_hull()
    logger.debug(shapely_shape)
    print(shapely_shape)


def build_concave_on_directory(input_dir):
    num_of_files_processed = 0
    num_of_files_failed = 0
    list_of_failed_files = []
    for filename in os.listdir(input_dir):

        extension = os.path.splitext(filename)[1]
        try:
            if(extension == ".nc"):
                filepath = os.path.join(input_dir, filename)
                ds = netCDF4.Dataset(filepath, 'r')
                ngu = NetCDFGridUtils(ds)
                shapely_shape = ngu.get_concave_hull()
                num_of_files_processed = num_of_files_processed + 1
                logger.debug("{} successful.".format(filename))
                del ds, ngu, shapely_shape

        except Exception as e:
            logger.debug("error on file: {}".format(filename))
            logger.debug(e)
            num_of_files_failed = num_of_files_failed + 1
            #list_of_failed_files.append(filename)

    logger.debug("Number of files proccessd: {}".format(num_of_files_processed))
    logger.debug("Number of files failed: {}".format(num_of_files_failed))

def main():

   # input_file = sys.argv[1]
    input_file = "C:/Users/u62231/Desktop/Projects/gadds/grid_examples/nonawagsgrids/ausbath_09_v4_ex_ex.nc"
   #input_file = "Bowen_Surat_Gravity_Bouguer_Anomaly.nc"
#GascoyneNorth_Complete_Sph_Cap_Bouguer_1VD_Geodetic.nc

    print(input_file)
    try:
        build_concave_on_single_file(input_file)
    except Exception as e:
        logger.debug(e)
        print(e)
    #build_concave_on_single_file()

if __name__ == "__main__":
    main()
