"""
Unit tests for geophys_utils._netcdf_grid_utils against a NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import os
import netCDF4
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils

netcdf_grid_utils = None

class TestNetCDFGridUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFGridUtils Constructor.
    N.B: This should be run first"""
    
    NC_PATH = 'test_grid.nc'
    
    def test_netcdf_grid_utils_constructor(self):
        print 'Testing NetCDFGridUtils constructor'
        global netcdf_grid_utils
        
        nc_path = os.path.join(os.path.dirname(__file__), TestNetCDFGridUtilsConstructor.NC_PATH)
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_grid_utils = NetCDFGridUtils(nc_dataset)
    
class TestNetCDFGridUtilsFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_grid_utils module."""
    
    MAX_BYTES = 1600
    MAX_ERROR = 0.000001
    TEST_COORDS = (148.213, -36.015)
    TEST_INDICES = [1, 1]
    TEST_FRACTIONAL_INDICES = [1.25, 1.25]
    TEST_VALUE = 0.0
    TEST_INTERPOLATED_VALUE = -99997.6171875
    
    def test_get_indices_from_coords(self):
        print 'Testing get_indices_from_coords function'
        indices = netcdf_grid_utils.get_indices_from_coords(TestNetCDFGridUtilsFunctions.TEST_COORDS)
        assert indices == TestNetCDFGridUtilsFunctions.TEST_INDICES, 'Indices incorrect'

    def test_get_fractional_indices_from_coords(self):
        print 'Testing get_fractional_indices_from_coords function'
        indices = netcdf_grid_utils.get_fractional_indices_from_coords(TestNetCDFGridUtilsFunctions.TEST_COORDS)
        for ordinate_index in range(len(indices)):
            assert round(indices[ordinate_index], 6) == TestNetCDFGridUtilsFunctions.TEST_FRACTIONAL_INDICES[ordinate_index], 'Fractional index incorrect'

    def test_get_value_at_coords(self):
        print 'Testing get_value_at_coords function'
        value = netcdf_grid_utils.get_value_at_coords(TestNetCDFGridUtilsFunctions.TEST_COORDS)
        assert value.data == TestNetCDFGridUtilsFunctions.TEST_VALUE, 'Incorrect retrieved value'

    def test_get_interpolated_value_at_coords(self):
        print 'Testing get_interpolated_value_at_coords function'
        interpolated_value = netcdf_grid_utils.get_interpolated_value_at_coords(TestNetCDFGridUtilsFunctions.TEST_COORDS)
        assert interpolated_value == TestNetCDFGridUtilsFunctions.TEST_INTERPOLATED_VALUE, 'Incorrect interpolated value retrieved'

    def test_sample_transect(self):
        print 'Testing sample_transect function'
        #TODO: Finish this!
        #transect_samples = netcdf_grid_utils.sample_transect(transect_vertices, crs=None, sample_metres=None):



# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFGridUtilsConstructor,
                    TestNetCDFGridUtilsFunctions]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())

if __name__ == '__main__':
    main()
