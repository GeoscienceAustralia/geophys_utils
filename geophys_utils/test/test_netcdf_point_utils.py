#!/usr/bin/env python

#===============================================================================
#    Copyright 2017 Geoscience Australia
# 
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#===============================================================================
"""
Unit tests for geophys_utils._netcdf_point_utils against a NetCDF line data file

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import os
import re
import netCDF4
import numpy as np
from geophys_utils._netcdf_point_utils import NetCDFPointUtils

netcdf_point_utils = None

#NC_PATH = '/g/data2/uc0/rr2_dev/axi547/GSSA_P1255MAG_Marree.nc'
NC_PATH = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/GSSA_P1255MAG_Marree.nc'
NC_TITLE = 'Marree Airborne Magnetic & Radiometric Survey, SA, 2012'
#NC_PATH = 'test_line.nc'
#NC_PATH = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/GSSA_P1255MAG_Marree.nc'

TEST_BOUNDS = (137, -29, 138, -28)
GRID_RESOLUTION = 0.001
RESAMPLING_METHOD = 'linear' # 'linear', 'nearest' or 'cubic'. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html

MAX_BYTES = 1600
MAX_ERROR = 0.000001
TEST_COORDS = (148.213, -36.015)
TEST_INDICES = [1, 1]
TEST_FRACTIONAL_INDICES = [1.25, 1.25]
TEST_VALUE = 0.0
TEST_INTERPOLATED_VALUE = -99997.6171875
SPATIAL_MASK_COUNT = 4613089


TEST_GRID_RESULTS = (('GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]]',
                      [136.4805, 0.001, 0, -27.988500000000002, 0, -0.001],
                      (1778, 3047)
                      ),
                     ('GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]]',
                      [136.9995, 0.001, 0, -27.9995, 0, -0.001],
                      (1001, 1001)
                      )
                     )

TEST_GET_LINE_MASK_RESULTS = ((190520, 5032),
                              (190500, 4994)
                              )

TEST_GET_LINE_RESULTS = ((190520, 4, 10064),
                         (190500, 4, 9988)
                         )

   
class TestNetCDFPointUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFPointUtils Constructor.
    N.B: This should be run first"""
    
    def test_netcdf_point_utils_constructor(self):
        print('Testing NetCDFPointUtils constructor')
        global netcdf_point_utils
        
        if re.match('^http.*', NC_PATH):
            nc_path = NC_PATH
        else:
            nc_path = os.path.join(os.path.dirname(__file__), NC_PATH)
        print(nc_path)   
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_point_utils = NetCDFPointUtils(nc_dataset)
        
        #print(netcdf_point_utils.__dict__)
        assert nc_dataset.title == NC_TITLE, 'Invalid dataset title: "{}" != "{}"'.format(nc_dataset.title, NC_TITLE)
    
class TestNetCDFPointUtilsFunctions1(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_point_utils functions"""
    
    def test_get_polygon(self):
        print('Testing get_polygon function')
        polygon = netcdf_point_utils.get_polygon()
        assert polygon is None, 'This is just plain messed up' #TODO: Find out why None is returned

    def test_get_spatial_mask(self):
        print('Testing get_spatial_mask function')
        spatial_mask = netcdf_point_utils.get_spatial_mask(TEST_BOUNDS)
        #print(spatial_mask)
        assert np.count_nonzero(spatial_mask) == SPATIAL_MASK_COUNT, 'Unexpected spatial mask count'
        

class TestNetCDFPointUtilsGridFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_point_utils functions"""
    
    def test_grid_points(self):
        print('Testing grid_points function')
        grids, crs, geotransform = netcdf_point_utils.grid_points(grid_resolution=GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 point_step = 100)
        assert (crs, geotransform, grids.shape) == TEST_GRID_RESULTS[0], 'Invalid grid results: {} != {}'.format((crs, geotransform, grids.shape), TEST_GRID_RESULTS[0])

        print('Testing bounded grid_points function')
        grids, crs, geotransform = netcdf_point_utils.grid_points(grid_resolution=GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 native_grid_bounds=TEST_BOUNDS,
                                                                 point_step = 100)
        assert (crs, geotransform, grids.shape) == TEST_GRID_RESULTS[1], 'Invalid grid results: {} != {}'.format((crs, geotransform, grids.shape), TEST_GRID_RESULTS[1])



# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFPointUtilsConstructor,
                    TestNetCDFPointUtilsFunctions1,
                    TestNetCDFPointUtilsGridFunctions
                    ]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())

if __name__ == '__main__':
    main()
