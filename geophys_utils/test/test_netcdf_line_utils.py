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

class TestNetCDFLineUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFLineUtils Constructor.
    N.B: This should be run first"""
    
    #NC_PATH = 'test_line.nc'
    NC_PATH = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/GSSA_P1255MAG_Marree.nc'
    
    def test_netcdf_line_utils_constructor(self):
        print 'Testing NetCDFLineUtils constructor'
        global netcdf_line_utils
        
        if re.match('^http.*', TestNetCDFLineUtilsConstructor.NC_PATH):
            nc_path = TestNetCDFLineUtilsConstructor.NC_PATH
        else:
            nc_path = os.path.join(os.path.dirname(__file__), TestNetCDFLineUtilsConstructor.NC_PATH)
        print nc_path    
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_line_utils = NetCDFLineUtils(nc_dataset)
        
        print netcdf_line_utils.__dict__
    
class TestNetCDFLineUtilsFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""
    
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
    
    def test_get_polygon(self):
        print 'Testing get_polygon function'
        polygon = netcdf_line_utils.get_polygon()
        print polygon

    def test_get_spatial_mask(self):
        print 'Testing get_spatial_mask function'
        spatial_mask = netcdf_line_utils.get_spatial_mask(TestNetCDFLineUtilsFunctions.TEST_BOUNDS)
        print spatial_mask
        print np.count_nonzero(spatial_mask)
        
    def test_get_line_masks(self):
        print 'Testing get_lines function'
        line_masks = netcdf_line_utils.get_line_masks()
        print line_masks
        
        for line_number in sorted(line_masks.keys()): 
            print line_number, np.count_nonzero(line_masks[line_number])
        
        line_masks = netcdf_line_utils.get_lines(bounds=TestNetCDFLineUtilsFunctions.TEST_BOUNDS)
        for line_number in sorted(line_masks.keys()): 
            print line_number, np.count_nonzero(line_masks[line_number])

class TestNetCDFLineUtilsGridFunctions(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""
    
    def test_grid_points(self):
        print 'Testing grid_points function'
        grids, crs, geotransform = netcdf_line_utils.grid_points(grid_resolution=TestNetCDFLineUtilsFunctions.GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 point_step = 100)
        print crs
        print geotransform
        print grids.shape

        grids, crs, geotransform = netcdf_line_utils.grid_points(grid_resolution=TestNetCDFLineUtilsFunctions.GRID_RESOLUTION, 
                                                                 variables='mag_awags',
                                                                 grid_bounds=TestNetCDFLineUtilsFunctions.TEST_BOUNDS,
                                                                 point_step = 100)
        print crs
        print geotransform
        print grids.shape



# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFLineUtilsConstructor,
                    TestNetCDFLineUtilsFunctions,
                    #TestNetCDFLineUtilsGridFunctions
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
