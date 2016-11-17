'''
Created on 16/11/2016

@author: Alex Ip
'''
import netCDF4
import numpy as np
import math
import os
import re
import tempfile
from scipy.interpolate import griddata
from geophys_utils._crs_utils import get_spatial_ref_from_crs, transform_coords

class NetCDFLineUtils(object):
    '''
    NetCDFLineUtils class to do various fiddly things with NetCDF geophysics line data files.
    '''

    def __init__(self, netcdf_dataset):
        '''
        NetCDFLineUtils Constructor
        '''
        self.netcdf_dataset = netcdf_dataset
        self.opendap = (re.match('^http.*', self.netcdf_dataset.filepath()) is not None)
        if self.opendap:
            self.max_bytes = 500000000 # 500MB limit for NCI OPeNDAP
        else:
            self.max_bytes = 8000000000 # 8GB limit for direct netCDF file access
        
        self.crs_variable = netcdf_dataset.variables['crs'] #TODO: Make this more general
        try:
            self.crs = self.crs_variable.spatial_ref
        except:
            self.crs = get_spatial_ref_from_crs(self.crs_variable.epsg_code).ExportToWkt()

        self.point_variables = list([var_name for var_name in self.netcdf_dataset.variables.keys() 
                                     if 'point' in self.netcdf_dataset.variables[var_name].dimensions
                                     and var_name not in ['latitude', 'longitude', 'point', 'fiducial', 'flag_linetype']
                                     ])
        
        # Create local cache for coordinates
        nc_cache_path = os.path.join(tempfile.gettempdir(), re.sub('\W', '_', os.path.splitext(self.netcdf_dataset.filepath())[0] + '.nc'))
        self._nc_cache_dataset = netCDF4.Dataset(nc_cache_path, mode="w", clobber=True, format='NETCDF4')
        
        point_dimension = self.netcdf_dataset.dimensions['point']
        self.point_count = len(point_dimension)
        
        self._nc_cache_dataset.createDimension('point', self.point_count if not point_dimension.isunlimited() else None)
        self._nc_cache_dataset.createDimension('yx', 2)
    
        var_options = {'zlib': False}
        self._nc_cache_dataset.createVariable('latlon', 
                                      self.netcdf_dataset.variables['longitude'].dtype, 
                                      ('point', 'yx'),
                                      **var_options
                                      )
        self.latlon = self._nc_cache_dataset.variables['latlon']

        pieces_required = max(self.latlon.dtype.itemsize * reduce(lambda x, y: x * y, self.latlon.shape) / self.max_bytes, 1)
        max_elements = self.point_count / pieces_required

        # Populate latlon cache
        start_index = 0
        while start_index < self.point_count:
            end_index = min(start_index + max_elements, self.point_count)
            # Pull lat/lon coordinates out of individual lat/lon arrays - note YX order
            latlon_slice = slice(start_index, end_index)
            self.latlon[latlon_slice, 0] = self.netcdf_dataset.variables['latitude'][latlon_slice]
            self.latlon[latlon_slice, 1] = self.netcdf_dataset.variables['longitude'][latlon_slice]
            start_index += max_elements
 
        # Determine exact spatial bounds
        min_lat = min(self.latlon[:,0])
        max_lat = max(self.latlon[:,0])
        min_lon = min(self.latlon[:,1])
        max_lon = max(self.latlon[:,1])
        self.bounds = (min_lon, min_lat, max_lon, max_lat)

        
    def __del__(self):
        cache_file_path = self._nc_cache_dataset.filename()
        self._nc_cache_dataset.close()
        os.remove(cache_file_path)
        
    def get_polygon(self):
        pass
    
    def get_spatial_mask(self, bounds):
        '''
        Return boolean mask of dimension 'point' for all coordinates within specified bounds
        '''
        return np.logical_and((bounds[1] <= self.latlon[:,0]) * (self.latlon[:,0] <= bounds[3]), 
                              (bounds[0] <= self.latlon[:,1]) * (self.latlon[:,1] <= bounds[2]))
    
    def get_lines(self, bounds=None, variables=None, lines=None):
        # Grid all data variables if not specified
        bounds = bounds or self.bounds
        variables = variables or self.point_variables
        lines = lines or self.netcdf_dataset.variables['line'][...]
    
    def grid_points(self, grid_resolution, grid_bounds=None, variables=None, resampling_method='linear', grid_crs=None):
        '''
        '''
        grid_bounds = grid_bounds or self.bounds

        # Allow single variable to be specified as a string not in a list
        if type(variables) in [str, unicode]:
            variables = [variables]

        # Extend area for points out beyond grid extents for nice interpolation at edges
        expanded_grid_bounds = [grid_bounds[0]-2*grid_resolution,
                                grid_bounds[1]-2*grid_resolution,
                                grid_bounds[2]+2*grid_resolution,
                                grid_bounds[3]+2*grid_resolution
                                ]

        spatial_subset_mask = self.get_spatial_mask(expanded_grid_bounds)
        
        # Transform grid extents
        if grid_crs is not None:
            grid_bounds = np.array(transform_coords(np.reshape(np.array(grid_bounds), (2.2)), self.crs, grid_crs)).flatten()

        # Grid all data variables if not specified
        variables = variables or self.point_variables
        
        # Determine spatial grid bounds rounded out to nearest GRID_RESOLUTION multiple
        min_lat = round(math.floor(grid_bounds[1] / grid_resolution) * grid_resolution, 6)
        max_lat = round(math.floor(grid_bounds[3] / grid_resolution + 1.0) * grid_resolution, 6)
        min_lon = round(math.floor(grid_bounds[0] / grid_resolution) * grid_resolution, 6)
        max_lon = round(math.floor(grid_bounds[2] / grid_resolution + 1.0) * grid_resolution, 6)
        grid_bounds = (min_lon, min_lat, max_lon, max_lat)    
        
    
        # Create grids of Y and X values. Note YX ordering and inverted Y
        # Note GRID_RESOLUTION/2 fudge to avoid truncation due to rounding error
        grid_y, grid_x = np.mgrid[grid_bounds[3]:grid_bounds[1]-grid_resolution/2:-grid_resolution, 
                                 grid_bounds[0]:grid_bounds[2]+grid_resolution/2:grid_resolution]

        # Skip points to reduce memory requirements - this is a horrible, lazy hack!
        #TODO: Implement function which grids spatial subsets.
        point_subset_mask = np.zeros(shape= self.netcdf_dataset.variables['point'].shape, dtype=bool)
        point_subset_mask[0:-1:2] = True
        point_subset_mask = np.logical_and(spatial_subset_mask, point_subset_mask)
        
        coordinates = self.latlon[...][point_subset_mask]
        # Reproject coordinates if required
        if grid_crs is not None:
            # N.B: Be careful about XY vs YX coordinate order         
            coordinates = transform_coords(coordinates[:,::-1], self.crs, grid_crs)[:,::-1]

        # Interpolate required values to the grid
        grids = {}
        for variable in [self.netcdf_dataset.variables[var_name] for var_name in variables]:
            grids[variable.name] = griddata(coordinates, 
                                  variable[...][point_subset_mask], #TODO: Check why this is faster than direct indexing
                                  (grid_y, grid_x), 
                                  method=resampling_method)
            
        return grids
