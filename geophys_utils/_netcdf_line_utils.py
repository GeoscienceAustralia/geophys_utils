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
from geophys_utils._crs_utils import get_spatial_ref_from_crs, transform_coords, get_utm_crs

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
        self._nc_cache_dataset.createDimension('xy', 2)
    
        var_options = {'zlib': False}
        self._nc_cache_dataset.createVariable('xycoords', 
                                      self.netcdf_dataset.variables['longitude'].dtype, 
                                      ('point', 'xy'),
                                      **var_options
                                      )
        self.xycoords = self._nc_cache_dataset.variables['xycoords']

        pieces_required = max(self.xycoords.dtype.itemsize * reduce(lambda x, y: x * y, self.xycoords.shape) / self.max_bytes, 1)
        max_elements = self.point_count / pieces_required

        # Populate xycoords cache
        start_index = 0
        while start_index < self.point_count:
            end_index = min(start_index + max_elements, self.point_count)
            # Pull lat/lon coordinates out of individual lat/lon arrays - note XY order
            xycoords_slice = slice(start_index, end_index)
            self.xycoords[xycoords_slice, 0] = self.netcdf_dataset.variables['longitude'][xycoords_slice]
            self.xycoords[xycoords_slice, 1] = self.netcdf_dataset.variables['latitude'][xycoords_slice]
            start_index += max_elements
 
        # Determine exact spatial bounds
        min_lon = min(self.xycoords[:,0])
        max_lon = max(self.xycoords[:,0])
        min_lat = min(self.xycoords[:,1])
        max_lat = max(self.xycoords[:,1])
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
        return np.logical_and((bounds[0] <= self.xycoords[:,0]) * (self.xycoords[:,0] <= bounds[2]), 
                              (bounds[1] <= self.xycoords[:,1]) * (self.xycoords[:,1] <= bounds[3]))
    
    def get_lines(self, bounds=None, variables=None, lines=None):
        bounds = bounds or self.bounds

        # Return all data variables if not specified
        variables = variables or self.point_variables
        
        # Allow single variable to be given as a string
        single_var = (type(variables) in [str, unicode])
        if single_var:
            variables = [variables]
        
        # Return all lines if not specified
        lines = lines or self.netcdf_dataset.variables['line'][...]
    
    def grid_points(self, grid_resolution, variables=None, native_grid_bounds=None, reprojected_grid_bounds=None, resampling_method='linear', grid_crs=None, point_step=1):
        '''
        '''
        def get_reprojected_bounds(bounds, from_crs, to_crs):
            if (to_crs is None) or (from_crs is None) or (to_crs == from_crs):
                return bounds
            
            original_bounding_box =((bounds[0], bounds[1]), (bounds[2], bounds[1]), (bounds[2], bounds[3]), (bounds[0], bounds[3]))
            reprojected_bounding_box = np.array(transform_coords(original_bounding_box, from_crs, to_crs))
            
            return [min(reprojected_bounding_box[:,0]), min(reprojected_bounding_box[:,1]), max(reprojected_bounding_box[:,0]), max(reprojected_bounding_box[:,1])]
            
            
        assert not (native_grid_bounds and reprojected_grid_bounds), 'Either native_grid_bounds or reprojected_grid_bounds can be provided, but not both'
        # Grid all data variables if not specified
        variables = variables or self.point_variables

        # Allow single variable to be given as a string
        single_var = (type(variables) in [str, unicode])
        if single_var:
            variables = [variables]
        
        if native_grid_bounds:
            reprojected_grid_bounds = get_reprojected_bounds(native_grid_bounds, self.crs, grid_crs)
        elif reprojected_grid_bounds:
            native_grid_bounds = get_reprojected_bounds(reprojected_grid_bounds, grid_crs, self.crs)
        else: # No reprojection required
            native_grid_bounds = self.bounds
            reprojected_grid_bounds = self.bounds

        print native_grid_bounds
        print reprojected_grid_bounds
        
        # Determine spatial grid bounds rounded out to nearest GRID_RESOLUTION multiple
        pixel_centre_bounds = (round(math.floor(reprojected_grid_bounds[0] / grid_resolution) * grid_resolution, 6),
                       round(math.floor(reprojected_grid_bounds[1] / grid_resolution) * grid_resolution, 6),
                       round(math.floor(reprojected_grid_bounds[2] / grid_resolution - 1.0) * grid_resolution + grid_resolution, 6),
                       round(math.floor(reprojected_grid_bounds[3] / grid_resolution - 1.0) * grid_resolution + grid_resolution, 6)
                       )
        
        grid_size = [pixel_centre_bounds[dim_index+2] - pixel_centre_bounds[dim_index] for dim_index in range(2)]

        # Extend area for points an arbitrary 4% out beyond grid extents for nice interpolation at edges
        expanded_grid_bounds = [pixel_centre_bounds[0]-grid_size[0]/50.0,
                                pixel_centre_bounds[1]-grid_size[0]/50.0,
                                pixel_centre_bounds[2]+grid_size[1]/50.0,
                                pixel_centre_bounds[3]+grid_size[1]/50.0
                                ]

        spatial_subset_mask = self.get_spatial_mask(get_reprojected_bounds(expanded_grid_bounds, grid_crs, self.crs))
        
        # Create grids of Y and X values. Note YX ordering and inverted Y
        # Note GRID_RESOLUTION/2.0 fudge to avoid truncation due to rounding error
        grid_y, grid_x = np.mgrid[pixel_centre_bounds[3]:pixel_centre_bounds[1]-grid_resolution/2.0:-grid_resolution, 
                                 pixel_centre_bounds[0]:pixel_centre_bounds[2]+grid_resolution/2.0:grid_resolution]

        # Skip points to reduce memory requirements
        #TODO: Implement function which grids spatial subsets.
        point_subset_mask = np.zeros(shape= self.netcdf_dataset.variables['point'].shape, dtype=bool)
        point_subset_mask[0:-1:point_step] = True
        point_subset_mask = np.logical_and(spatial_subset_mask, point_subset_mask)
        
        coordinates = self.xycoords[...][point_subset_mask]
        # Reproject coordinates if required
        if grid_crs is not None:
            # N.B: Be careful about XY vs YX coordinate order         
            coordinates = np.array(transform_coords(coordinates[...], self.crs, grid_crs))

        # Interpolate required values to the grid - Note YX ordering for image
        grids = {}
        for variable in [self.netcdf_dataset.variables[var_name] for var_name in variables]:
            grids[variable.name] = griddata(coordinates[:,::-1],
                                  variable[...][point_subset_mask], #TODO: Check why this is faster than direct indexing
                                  (grid_y, grid_x), 
                                  method=resampling_method)

        if single_var:
            grids = grids.values()[0]
            
        #  crs:GeoTransform = "109.1002342895272 0.00833333 0 -9.354948067227777 0 -0.00833333 "
        geotransform = [pixel_centre_bounds[0]-grid_resolution/2.0,
                        grid_resolution,
                        0,
                        pixel_centre_bounds[3]+grid_resolution/2.0,
                        0,
                        -grid_resolution
                        ] 

        return (grid_crs or self.crs), geotransform, grids
    
    
    def utm_grid_points(self, utm_grid_resolution, variables=None, native_grid_bounds=None, resampling_method='linear', point_step=1):
        native_grid_bounds = native_grid_bounds or self.bounds
        
        native_centre_coords = [(native_grid_bounds[dim_index] + native_grid_bounds[dim_index+2]) / 2.0 for dim_index in range(2)]
        utm_crs = get_utm_crs(native_centre_coords, self.crs)
        
        return self.grid_points(grid_resolution=utm_grid_resolution, 
                                variables=variables,
                                native_grid_bounds=native_grid_bounds, 
                                resampling_method=resampling_method, 
                                grid_crs=utm_crs, 
                                point_step=point_step
                                )

