'''
Created on 16/11/2016

@author: Alex Ip
'''
import numpy as np
import math
from scipy.interpolate import griddata
from geophys_utils._crs_utils import get_spatial_ref_from_crs

class NetCDFLineUtils(object):
    '''
    NetCDFLineUtils class to do various fiddly things with NetCDF geophysics line data files.
    '''

    def __init__(self, netcdf_dataset):
        '''
        NetCDFLineUtils Constructor
        '''
        self.netcdf_dataset = netcdf_dataset
        
        self.crs_variable = netcdf_dataset.variables['crs'] #TODO: Make this more general
        try:
            self.crs = self.crs_variable.spatial_ref
        except:
            get_spatial_ref_from_crs(self.crs_variable.epsg_code).ExportToWkt()

        
        # Pull lat/lon coordinates out of individual lat/lon arrays - note YX order
        self.coord_array = np.zeros(shape=(self.netcdf_dataset.variables['point'].shape[0], 2), 
                               dtype=self.netcdf_dataset.variables['longitude'].dtype)
        self.coord_array[:,0] = self.netcdf_dataset.variables['latitude'][:]
        self.coord_array[:,1] = self.netcdf_dataset.variables['longitude'][:]
 
        # Determine exact spatial bounds
        min_lat = min(self.coord_array[:,0])
        max_lat = max(self.coord_array[:,0])
        min_lon = min(self.coord_array[:,1])
        max_lon = max(self.coord_array[:,1])
        self.bounds = (min_lon, min_lat, max_lon, max_lat)

        
    def get_polygon(self):
        pass
    
    def get_spatial_mask(self, bounds):
        '''
        Return boolean mask of dimension 'point' for all coordinates within specified bounds
        '''
        return (bounds[0] <= self.coord_array[:,0] <= bounds[2]) and (bounds[1] <= self.coord_array[:,1] <= bounds[3])
    
    def get_lines(self, bounds, variables):
        pass
    
    def grid_points(self, grid_bounds, grid_resolution, variables=None, resampling_method='linear', grid_crs=None):
        '''
        '''
        grid_bounds = grid_bounds or self.bounds
        expanded_grid_bounds = [grid_bounds[0]-2*grid_resolution,
                                grid_bounds[1]-2*grid_resolution,
                                grid_bounds[2]+2*grid_resolution,
                                grid_bounds[3]+2*grid_resolution
                                ]
        spatial_subset_mask = self.get_spatial_mask(expanded_grid_bounds)
        
        # Grid all data variables if not specified
        variables = variables or [var_name for var_name in self.netcdf_dataset.variables.keys() 
                                  if 'point' in self.netcdf_dataset.variables[var_name].dimensions
                                  and var_name not in ['latitude', 'longitude', 'point', 'fiducial', 'flag_linetype']
                                  ]
        
        # Determine spatial grid bounds rounded out to nearest GRID_RESOLUTION multiple
        min_lat = round(math.floor(grid_bounds[0] / grid_resolution) * grid_resolution, 6)
        max_lat = round(math.floor(grid_bounds[0] / grid_resolution + 1.0) * grid_resolution, 6)
        min_lon = round(math.floor(grid_bounds[0] / grid_resolution) * grid_resolution, 6)
        max_lon = round(math.floor(grid_bounds[0] / grid_resolution + 1.0) * grid_resolution, 6)
        grid_bounds = (min_lon, min_lat, max_lon, max_lat)            
    
        # Create grids of Y and X values. Note YX ordering and inverted Y
        # Note GRID_RESOLUTION/2 fudge to avoid truncation due to rounding error
        grid_y, grid_x = np.mgrid[grid_bounds[3]:grid_bounds[1]-grid_resolution/2:-grid_resolution, 
                                 grid_bounds[0]:grid_bounds[2]+grid_resolution/2:grid_resolution]

        # Skip points to reduce memory requirements - this is a horrible, lazy hack!
        #TODO: Implement function which grids spatial subsets.
        point_subset_mask = np.zeros(shape= self.netcdf_dataset.variables['point'].shape, dtype=bool)
        point_subset_mask[0:-1:2] = True
        point_subset_mask = spatial_subset_mask and point_subset_mask
        
        # Interpolate required values to the grid
        grids = {}
        for variable in [self.netcdf_dataset.variables[var_name] for var_name in variables]:
            grids[variable.name] = griddata(self.coord_array[point_subset_mask], 
                                  variable[point_subset_mask], 
                                  (grid_y, grid_x), 
                                  method=resampling_method)
            
        return grids

    def clean_lines(self, lines):
        pass
