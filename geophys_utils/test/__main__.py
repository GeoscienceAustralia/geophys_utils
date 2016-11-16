"""
Main unit for test module
Unit tests for ncskosdump and ld_functions against a modified NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
from geophys_utils.test import test_array_pieces, test_crs_utils, test_data_stats, test_netcdf_grid_utils

# Run all tests
test_array_pieces.main()
test_crs_utils.main()
test_data_stats.main()
test_netcdf_grid_utils.main()