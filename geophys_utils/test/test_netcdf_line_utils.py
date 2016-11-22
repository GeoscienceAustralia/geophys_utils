"""
Unit tests for geophys_utils._netcdf_line_utils against a NetCDF line data file

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import os
import re
import netCDF4
import numpy as np
from geophys_utils._netcdf_line_utils import NetCDFLineUtils

netcdf_line_utils = None

#NC_PATH = '/g/data2/uc0/rr2_dev/axi547/GSSA_P1255MAG_Marree.nc'
NC_PATH = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/GSSA_P1255MAG_Marree.nc'
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
    
class TestNetCDFLineUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFLineUtils Constructor.
    N.B: This should be run first"""
    
    def test_netcdf_line_utils_constructor(self):
        print 'Testing NetCDFLineUtils constructor'
        global netcdf_line_utils
        
        if re.match('^http.*', NC_PATH):
            nc_path = NC_PATH
        else:
            nc_path = os.path.join(os.path.dirname(__file__), NC_PATH)
        print nc_path    
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_line_utils = NetCDFLineUtils(nc_dataset)
        
        print netcdf_line_utils.__dict__
    
class TestNetCDFLineUtilsFunctions1(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""
    
    def test_get_polygon(self):
        print 'Testing get_polygon function'
        polygon = netcdf_line_utils.get_polygon()
        print polygon

    def test_get_spatial_mask(self):
        print 'Testing get_spatial_mask function'
        spatial_mask = netcdf_line_utils.get_spatial_mask(TEST_BOUNDS)
        print spatial_mask
        print np.count_nonzero(spatial_mask)
        
    def test_get_line_masks(self):
        print 'Testing get_lines function'
        for line_number, line_mask in netcdf_line_utils.get_line_masks():
            print 'Line %d has %d points' % (line_number, np.count_nonzero(line_mask))


class TestNetCDFLineUtilsFunctions2(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""
    
    def test_get_lines(self):
        print 'Testing get_lines function'
        for line_number, line_dict in netcdf_line_utils.get_lines():
            print 'Line %d has %d variables with %d points' % (line_number,
                                                               len(line_dict)-1, 
                                                               np.count_nonzero(line_dict['coordinates'])
                                                               )


class TestNetCDFLineUtilsGridFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""
    
    def test_grid_points(self):
        print 'Testing grid_points function'
        grids, crs, geotransform = netcdf_line_utils.grid_points(grid_resolution=GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 point_step = 100)
        print crs
        print geotransform
        print grids.shape

        grids, crs, geotransform = netcdf_line_utils.grid_points(grid_resolution=GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 native_grid_bounds=TEST_BOUNDS,
                                                                 point_step = 100)
        print crs
        print geotransform
        print grids.shape



# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFLineUtilsConstructor,
                    TestNetCDFLineUtilsFunctions1,
                    TestNetCDFLineUtilsFunctions2,
                    TestNetCDFLineUtilsGridFunctions
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
