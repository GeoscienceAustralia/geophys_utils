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
'''
Created on 16Nov.,2016

@author: u76345
'''
import logging
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

#===============================================================================
# if not logger.handlers:
#     # Set handler for root logger to standard output
#     console_handler = logging.StreamHandler(sys.stdout)
#     #console_handler.setLevel(logging.INFO)
#     console_handler.setLevel(logging.DEBUG)
#     console_formatter = logging.Formatter('%(message)s')
#     console_handler.setFormatter(console_formatter)
#     logger.addHandler(console_handler)
#===============================================================================
        
from geophys_utils._netcdf_utils import NetCDFUtils
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils
from geophys_utils._netcdf_point_utils import NetCDFPointUtils
from geophys_utils._netcdf_line_utils import NetCDFLineUtils
from geophys_utils._csw_utils import CSWUtils
from geophys_utils._array_pieces import array_pieces
from geophys_utils._data_stats import DataStats
from geophys_utils._polygon_utils import get_grid_edge_points, get_netcdf_edge_points, points2convex_hull, points2alpha_shape, netcdf2convex_hull
from geophys_utils._crs_utils import get_spatial_ref_from_wkt, get_coordinate_transformation, get_utm_wkt, transform_coords
from geophys_utils._transect_utils import line_length, point_along_line, utm_coords, coords2distance, sample_transect
from geophys_utils._dem_utils import DEMUtils
from geophys_utils._datetime_utils import date_string2datetime

try:
    from geophys_utils._gdal_grid_utils import get_gdal_wcs_dataset, get_gdal_grid_values
    from geophys_utils._array2file import array2file
except:
    logger.warning('Unable to import get_gdal_wcs_dataset, get_gdal_grid_values or array2file (no GDAL available')
