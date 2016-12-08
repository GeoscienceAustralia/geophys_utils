'''
Created on 16Nov.,2016

@author: u76345
'''
from geophys_utils._array_pieces import array_pieces
from geophys_utils._data_stats import DataStats
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils
from geophys_utils._netcdf_line_utils import NetCDFLineUtils
from geophys_utils._polygon_utils import get_grid_edge_points, get_netcdf_edge_points, points2convex_hull, points2alpha_shape, netcdf2convex_hull
from geophys_utils._crs_utils import get_spatial_ref_from_crs, get_coordinate_transformation, get_utm_crs, transform_coords
from geophys_utils._gdal_grid_utils import get_gdal_wcs_dataset, get_gdal_grid_values