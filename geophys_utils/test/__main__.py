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
Main unit for test module
Unit tests for ncskosdump and ld_functions against a modified NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
from geophys_utils.test import \
    test_array_pieces, \
    test_crs_utils, \
    test_data_stats, \
    test_netcdf_grid_utils, \
    test_netcdf_line_utils, \
    test_netcdf_point_utils

# Run all tests
test_array_pieces.main()
test_crs_utils.main()
test_data_stats.main()
test_netcdf_grid_utils.main()
test_netcdf_line_utils.main()
test_netcdf_point_utils.main()
