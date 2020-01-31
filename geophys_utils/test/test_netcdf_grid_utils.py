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
Unit tests for geophys_utils._netcdf_grid_utils against a NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import os
import netCDF4
import numpy as np
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils
from shapely.geometry.polygon import Polygon

netcdf_grid_utils = None

NC_PATH = 'test_grid.nc'    
MAX_BYTES = 1600
MAX_ERROR = 0.000001
TEST_COORDS = (148.213, -36.015)
TEST_MULTI_COORDS = np.array([[148.213, -36.015], [148.516, -35.316]])
TEST_INDICES = [1, 1]
TEST_MULTI_INDICES = [[1, 1], [176, 77]]
TEST_FRACTIONAL_INDICES = [1.25, 1.25]
TEST_VALUE = -99999.
TEST_MULTI_VALUES = [-99999.0, -134.711334229]
TEST_INTERPOLATED_VALUE = -99997.6171875
    
class TestNetCDFGridUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFGridUtils Constructor.
    N.B: This should be run first"""
    
    def test_netcdf_grid_utils_constructor(self):
        print('Testing NetCDFGridUtils constructor')
        global netcdf_grid_utils
        
        nc_path = os.path.join(os.path.dirname(__file__), NC_PATH)
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_grid_utils = NetCDFGridUtils(nc_dataset)
    
class TestNetCDFGridUtilsFunctions1(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_grid_utils functions"""
    
    def test_get_indices_from_coords(self):
        print('Testing get_indices_from_coords function with single coordinate {}'.format(TEST_COORDS))
        indices = netcdf_grid_utils.get_indices_from_coords(TEST_COORDS)
        assert (indices == np.array(TEST_INDICES)).all, 'Incorrect indices: {} instead of {}'.format(indices, TEST_INDICES)

        print('Testing get_indices_from_coords function with multi coordinates {}'.format(TEST_MULTI_COORDS))
        multi_indices = netcdf_grid_utils.get_indices_from_coords(TEST_MULTI_COORDS)
        assert (multi_indices == np.array(TEST_MULTI_INDICES)).all, 'Incorrect indices: {} instead of {}'.format(multi_indices, TEST_MULTI_INDICES)

    def test_get_fractional_indices_from_coords(self):
        print('Testing get_fractional_indices_from_coords function')
        indices = netcdf_grid_utils.get_fractional_indices_from_coords(TEST_COORDS)
        for ordinate_index in range(len(indices)):
            assert round(indices[ordinate_index], 6) == TEST_FRACTIONAL_INDICES[ordinate_index], 'Fractional index incorrect'

    def test_get_value_at_coords(self):
        print('Testing get_value_at_coords function with single coordinate {}'.format(TEST_COORDS))
        value = netcdf_grid_utils.get_value_at_coords(TEST_COORDS)
        assert value[0] == TEST_VALUE, 'Incorrect retrieved value: {} instead of {}'.format(value.data, TEST_VALUE)
        
        print('Testing get_value_at_coords function with multiple coordinates {}'.format(TEST_MULTI_COORDS))
        multi_values = netcdf_grid_utils.get_value_at_coords(TEST_MULTI_COORDS)
        assert (np.abs(np.array(multi_values) - np.array(TEST_MULTI_VALUES)) < MAX_ERROR).all(), 'Incorrect retrieved value: {} instead of {}'.format(multi_values, TEST_MULTI_VALUES)

    def test_get_interpolated_value_at_coords(self):
        print('Testing get_interpolated_value_at_coords function')
        interpolated_value = netcdf_grid_utils.get_interpolated_value_at_coords(TEST_COORDS)
        assert interpolated_value == TEST_INTERPOLATED_VALUE, 'Incorrect interpolated value retrieved'

    def test_sample_transect(self):
        print('Testing sample_transect function')
        #TODO: Finish this!
        #transect_samples = netcdf_grid_utils.sample_transect(transect_vertices, crs=None, sample_metres=None)

    def test_concave_hull(self):
        print('Testing concave hull')
        assert isinstance(netcdf_grid_utils.get_concave_hull(), Polygon)


# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFGridUtilsConstructor,
                    TestNetCDFGridUtilsFunctions1]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())

if __name__ == '__main__':
    main()
