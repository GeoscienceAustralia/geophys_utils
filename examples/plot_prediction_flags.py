import netCDF4
from geophys_utils import NetCDFPointUtils
from plot_point_dataset import plot_point_dataset
from pprint import pprint
import re
import numpy as np
from geophys_utils._transect_utils import utm_coords, coords2distance
import os
import sys
from geophys_utils import CSWUtils
from geophys_utils import DataStats # https://github.com/alex-ip/geophys2netcdf/tree/develop
from netCDF4 import Dataset
from pprint import pprint
import sqlite3


def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)

    return None


#netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/aem/AusAEM_Year1_Tranche1_Final_EM.nc'
#netcdf_path = "D:\\Temp\\AEM Data\\AusAEM_Year1_Tranche1_Final_EM.nc"
#netcdf_path = "C:\\Users\\u62231\\Desktop\\bad_line_datasets\\new\\GSNSW_P747MAG.nc"
#C:\\Users\\u62231\\Desktop\\bad_line_datasets\\new\\
input_netcdf_folder_path = sys.argv[1]

database = "C:\\Users\\u62231\\Desktop\\sqlite\\dataset_metadata_cache.sqlite"

conn = create_connection(database)
cur = conn.cursor()
rows = cur.fetchall()

for row in rows:
    print(row)

def main():

    for filename in os.listdir(input_netcdf_folder_path):
        if filename.endswith(".nc"):
            #if re.search('P861_862_870_871MAG', filename):
                print(filename)
                netcdf_dataset = netCDF4.Dataset("{}\{}".format(input_netcdf_folder_path, filename))
                npu = NetCDFPointUtils(netcdf_dataset)


                coord_predict = netcdf_dataset['coord_predicted'][:]
                coord_predict2 = np.ma.filled(coord_predict.astype(float), 9)

                unique, counts = np.unique(coord_predict2, return_counts=True)
                dictionary = dict(zip(unique, counts))

                filename = re.sub('.nc', '', filename)
                if 1 in dictionary or 2 in dictionary: # dictionary[1] is the key for interpolated values, 2 is extrapolated
                    try:
                        # Plot subset of survey
                        plot_point_dataset(npu,
                                           'coord_predicted',
                                           plot_title=filename,
                                           colour_scheme='jet',
                                           save_path='C:\\Users\\u62231\\Desktop\\predicted_value_figures\\' + filename + ".png",
                                           count_dict=dictionary)
                                           #save_path="C:\\Users\\u62231\\Desktop\\predicted_value_figures\\GeorginaBasinNTQLD1964.png")
                    except:
                        pass
