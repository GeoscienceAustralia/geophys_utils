#!/usr/bin/env python

# ===============================================================================
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
# ===============================================================================
"""
Unit tests for geophys_utils._netcdf_point_utils against a NetCDF line data file

Created on 15/11/2016

@author: Alex Ip
"""
import os
import re
import unittest

import netCDF4
import numpy as np
from shapely.geometry.polygon import Polygon

from geophys_utils._netcdf_point_utils import NetCDFPointUtils

netcdf_point_utils = None

NC_PATH = 'https://dapds00.nci.org.au/thredds/dodsC/iv65/Geoscience_Australia_Geophysics_Reference_Data_Collection/ground_gravity/SA/point/P194830/P194830-point-gravity.nc'
NC_TITLE = 'Leigh Creek Gravity (P194830), gravity point data'
TEST_BOUNDS = (138.23, -30.3, 138.45, -30.125)
GRID_RESOLUTION = 0.001
RESAMPLING_METHOD = 'linear'  # 'linear', 'nearest' or 'cubic'. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
VARIABLES = 'bouguer'
POINT_STEP = 1

MAX_BYTES = 1600
MAX_ERROR = 0.000001
TEST_COORDS = (148.213, -36.015)
TEST_INDICES = [1, 1]
TEST_FRACTIONAL_INDICES = [1.25, 1.25]
TEST_VALUE = 0.0
TEST_INTERPOLATED_VALUE = -99997.6171875
SPATIAL_MASK_COUNT = 158

TEST_GRID_RESULTS = (
    (
        'GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]]',
        [138.29250000000002, 0.001, 0, -30.1185, 0, -0.001],
        (172, 168)
    ),
    (
        'GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]]',
        [138.2295, 0.001, 0, -30.1245, 0, -0.001],
        (176, 221)
    )
)

EXPECTED_GML = 'POLYGON((138.3498 -30.2894, 138.3093 -30.1768, 138.2940 -30.1189, 138.3097 -30.1187, 138.3189 -30.1187, 138.3684 -30.1188, 138.4608 -30.2419, 138.4319 -30.2558, 138.3663 -30.2854, 138.3582 -30.2875, 138.3498 -30.2894))'


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

        # print(netcdf_point_utils.__dict__)
        assert nc_dataset.title == NC_TITLE, 'Invalid dataset title: "{}" != "{}"'.format(nc_dataset.title, NC_TITLE)


class TestNetCDFPointUtilsFunctions1(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_point_utils functions"""

    def test_get_polygon(self):
        print('Testing get_polygon function')
        polygon = netcdf_point_utils.get_polygon()
        assert polygon == EXPECTED_GML

    def test_get_spatial_mask(self):
        print('Testing get_spatial_mask function')
        spatial_mask = netcdf_point_utils.get_spatial_mask(TEST_BOUNDS)
        # print(spatial_mask)
        assert np.count_nonzero(
            spatial_mask) == SPATIAL_MASK_COUNT, f'Unexpected spatial mask count {np.count_nonzero(spatial_mask)} != {SPATIAL_MASK_COUNT}'

    def test_concave_hull(self):
        print('Testing concave hull')
        assert isinstance(netcdf_point_utils.get_concave_hull(), Polygon)


class TestNetCDFPointUtilsGridFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_point_utils functions"""

    def test_grid_points(self):
        print('Testing grid_points function')
        grids, crs, geotransform = netcdf_point_utils.grid_points(grid_resolution=GRID_RESOLUTION,
                                                                  variables=VARIABLES,
                                                                  point_step=POINT_STEP)

        assert (crs, geotransform, grids.shape) == TEST_GRID_RESULTS[0], 'Invalid grid results: {} != {}'.format(
            (crs, geotransform, grids.shape), TEST_GRID_RESULTS[0])

    def test_grid_points_bounded(self):
        print('Testing bounded grid_points function')
        grids, crs, geotransform = netcdf_point_utils.grid_points(grid_resolution=GRID_RESOLUTION,
                                                                  variables=VARIABLES,
                                                                  native_grid_bounds=TEST_BOUNDS,
                                                                  point_step=POINT_STEP)
        assert (crs, geotransform, grids.shape) == TEST_GRID_RESULTS[1], 'Invalid grid results: {} != {}'.format(
            (crs, geotransform, grids.shape), TEST_GRID_RESULTS[1])


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
