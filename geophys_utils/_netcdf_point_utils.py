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
from geophys_utils._crs_utils import transform_coords, get_utm_wkt
from geophys_utils._transect_utils import utm_coords, coords2distance
from geophys_utils._netcdf_utils import NetCDFUtils
from scipy.spatial.ckdtree import cKDTree
import logging

# Setup logging handlers if required
logger = logging.getLogger(__name__) # Get __main__ logger
logger.setLevel(logging.INFO) # Initial logging level for this module
    
class NetCDFPointUtils(NetCDFUtils):
    '''
    NetCDFPointUtils class to do various fiddly things with NetCDF geophysics point data files.
    '''

    def __init__(self, netcdf_dataset, debug=False):
        '''
        NetCDFPointUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a line dataset
        '''
        # Start of init function - Call inherited constructor first
        NetCDFUtils.__init__(self, netcdf_dataset, debug)

        self.point_variables = list([var_name for var_name in self.netcdf_dataset.variables.keys() 
                                     if 'point' in self.netcdf_dataset.variables[var_name].dimensions
                                     and var_name not in ['latitude', 'longitude', 'point', 'fiducial', 'flag_linetype']
                                     ])
        
        # Create local cache for coordinates
        nc_cache_path = os.path.join(tempfile.gettempdir(), re.sub('\W', '_', os.path.splitext(self.netcdf_dataset.filepath())[0] + '.nc'))
        self._nc_cache_dataset = netCDF4.Dataset(nc_cache_path, mode="w", clobber=True, format='NETCDF4')
        
        point_dimension = self.netcdf_dataset.dimensions['point']
        self.point_count = len(point_dimension)
        self.unlimited_points = point_dimension.isunlimited()
        
        self._nc_cache_dataset.createDimension('point', self.point_count if not self.unlimited_points else None)
        self._nc_cache_dataset.createDimension('xy', 2)
        
        var_options = self.netcdf_dataset.variables['longitude'].filters() or {}
        var_options['zlib'] = False
        if hasattr(self.netcdf_dataset.variables['longitude'], '_FillValue'):
            var_options['fill_value'] = self.netcdf_dataset.variables['longitude']._FillValue

        self._nc_cache_dataset.createVariable('xycoords', 
                                      self.netcdf_dataset.variables['longitude'].dtype, 
                                      ('point', 'xy'),
                                      **var_options
                                      )
        self.xycoords = self._nc_cache_dataset.variables['xycoords']
        self.xycoords[:,0] = self.fetch_array(self.netcdf_dataset.variables['longitude'])
        self.xycoords[:,1] = self.fetch_array(self.netcdf_dataset.variables['latitude'])
 
        # Determine exact spatial bounds
        min_lon = np.nanmin(self.xycoords[:,0])
        max_lon = np.nanmax(self.xycoords[:,0])
        min_lat = np.nanmin(self.xycoords[:,1])
        max_lat = np.nanmax(self.xycoords[:,1])
        
        # Create nested list of bounding box corner coordinates
        self.native_bbox = [[min_lon, min_lat], [max_lon, min_lat], [max_lon, max_lat], [min_lon, max_lat]]
        self.wgs84_bbox = transform_coords(self.native_bbox, from_wkt=self.wkt, to_wkt='EPSG:4326')

        # Define bounds
        self.bounds = [min_lon, min_lat, max_lon, max_lat]
               
        self.kdtree = None
        
    def __del__(self):
        '''
        NetCDFPointUtils Destructor
        '''
        try:
            cache_file_path = self._nc_cache_dataset.filepath()
            self._nc_cache_dataset.close()
            os.remove(cache_file_path)
        except:
            pass
        
    def fetch_array(self, source_array):
        '''
        Helper function to retrieve entire 1D array in pieces < self.max_bytes in size
        '''
        source_len = source_array.shape[0]
        pieces_required = max(source_array[0].itemsize * source_len // self.max_bytes, 1)
        max_elements = source_len // pieces_required
        
        cache_array = np.zeros((source_len,), dtype=source_array.dtype)

        # Copy array in pieces
        start_index = 0
        while start_index < source_len:
            end_index = min(start_index + max_elements, source_len)
            array_slice = slice(start_index, end_index)
            cache_array[array_slice] = source_array[array_slice]
            cache_array[array_slice] = source_array[array_slice]
            start_index += max_elements
            
        return cache_array
        
    def get_polygon(self):
        '''
        Under construction - do not use except for testing
        '''
        pass
    
    def get_spatial_mask(self, bounds, bounds_wkt=None):
        '''
        Return boolean mask of dimension 'point' for all coordinates within specified bounds and CRS
        '''
        if bounds_wkt is None:
            coordinates = self.xycoords
        else:
            coordinates = np.array(transform_coords(self.xycoords[...], self.wkt, bounds_wkt))
            
        return np.logical_and(np.logical_and((bounds[0] <= coordinates[:,0]), (coordinates[:,0] <= bounds[2])), 
                              np.logical_and((bounds[1] <= coordinates[:,1]), (coordinates[:,1] <= bounds[3]))
                              )
            
        
    
    def get_reprojected_bounds(self, bounds, from_wkt, to_wkt):
        '''
        Function to take a bounding box specified in one CRS and return its smallest containing bounding box in a new CRS
        @parameter bounds: bounding box specified as tuple(xmin, ymin, xmax, ymax) in CRS from_wkt
        @parameter from_wkt: WKT for CRS from which to transform bounds
        @parameter to_wkt: WKT for CRS to which to transform bounds
        
        @return reprojected_bounding_box: bounding box specified as tuple(xmin, ymin, xmax, ymax) in CRS to_wkt
        '''
        if (to_wkt is None) or (from_wkt is None) or (to_wkt == from_wkt):
            return bounds
        
        # Need to look at all four bounding box corners, not just LL & UR
        original_bounding_box =((bounds[0], bounds[1]), (bounds[2], bounds[1]), (bounds[2], bounds[3]), (bounds[0], bounds[3]))
        reprojected_bounding_box = np.array(transform_coords(original_bounding_box, from_wkt, to_wkt))
        
        return [min(reprojected_bounding_box[:,0]), min(reprojected_bounding_box[:,1]), max(reprojected_bounding_box[:,0]), max(reprojected_bounding_box[:,1])]
            
            
    def grid_points(self, grid_resolution, 
                    variables=None, 
                    native_grid_bounds=None, 
                    reprojected_grid_bounds=None, 
                    resampling_method='linear', 
                    grid_wkt=None, 
                    point_step=1):
        '''
        Function to grid points in a specified bounding rectangle to a regular grid of the specified resolution and crs
        @parameter grid_resolution: cell size of regular grid in grid CRS units
        @parameter variables: Single variable name string or list of multiple variable name strings. Defaults to all point variables
        @parameter native_grid_bounds: Spatial bounding box of area to grid in native coordinates 
        @parameter reprojected_grid_bounds: Spatial bounding box of area to grid in grid coordinates
        @parameter resampling_method: Resampling method for gridding. 'linear' (default), 'nearest' or 'cubic'. 
        See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html 
        @parameter grid_wkt: WKT for grid coordinate reference system. Defaults to native CRS
        @parameter point_step: Sampling spacing for points. 1 (default) means every point, 2 means every second point, etc.
        
        @return grids: dict of grid arrays keyed by variable name if parameter 'variables' value was a list, or
        a single grid array if 'variable' parameter value was a string
        @return wkt: WKT for grid coordinate reference system.
        @return geotransform: GDAL GeoTransform for grid
        '''
        assert not (native_grid_bounds and reprojected_grid_bounds), 'Either native_grid_bounds or reprojected_grid_bounds can be provided, but not both'
        # Grid all data variables if not specified
        variables = variables or self.point_variables

        # Allow single variable to be given as a string
        single_var = (type(variables) == str)
        if single_var:
            variables = [variables]
        
        if native_grid_bounds:
            reprojected_grid_bounds = self.get_reprojected_bounds(native_grid_bounds, self.wkt, grid_wkt)
        elif reprojected_grid_bounds:
            native_grid_bounds = self.get_reprojected_bounds(reprojected_grid_bounds, grid_wkt, self.wkt)
        else: # No reprojection required
            native_grid_bounds = self.bounds
            reprojected_grid_bounds = self.bounds

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

        spatial_subset_mask = self.get_spatial_mask(self.get_reprojected_bounds(expanded_grid_bounds, grid_wkt, self.wkt))
        
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
        if grid_wkt is not None:
            # N.B: Be careful about XY vs YX coordinate order         
            coordinates = np.array(transform_coords(coordinates[...], self.wkt, grid_wkt))

        # Interpolate required values to the grid - Note YX ordering for image
        grids = {}
        for variable in [self.netcdf_dataset.variables[var_name] for var_name in variables]:
            grids[variable.name] = griddata(coordinates[:,::-1],
                                  variable[...][point_subset_mask], #TODO: Check why this is faster than direct indexing
                                  (grid_y, grid_x), 
                                  method=resampling_method)

        if single_var:
            grids = list(grids.values())[0]
            
        #  crs:GeoTransform = "109.1002342895272 0.00833333 0 -9.354948067227777 0 -0.00833333 "
        geotransform = [pixel_centre_bounds[0]-grid_resolution/2.0,
                        grid_resolution,
                        0,
                        pixel_centre_bounds[3]+grid_resolution/2.0,
                        0,
                        -grid_resolution
                        ] 

        return grids, (grid_wkt or self.wkt), geotransform
    
    
    def utm_grid_points(self, utm_grid_resolution, variables=None, native_grid_bounds=None, resampling_method='linear', point_step=1):
        '''
        Function to grid points in a specified native bounding rectangle to a regular grid of the specified resolution in its local UTM CRS
        @parameter grid_resolution: cell size of regular grid in metres (UTM units)
        @parameter variables: Single variable name string or list of multiple variable name strings. Defaults to all point variables
        @parameter native_grid_bounds: Spatial bounding box of area to grid in native coordinates 
        @parameter resampling_method: Resampling method for gridding. 'linear' (default), 'nearest' or 'cubic'. 
        See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html 
        @parameter grid_wkt: WKT for grid coordinate reference system. Defaults to native CRS
        @parameter point_step: Sampling spacing for points. 1 (default) means every point, 2 means every second point, etc.
        
        @return grids: dict of grid arrays keyed by variable name if parameter 'variables' value was a list, or
        a single grid array if 'variable' parameter value was a string
        @return wkt: WKT for grid coordinate reference system (i.e. local UTM zone)
        @return geotransform: GDAL GeoTransform for grid
        '''
        native_grid_bounds = native_grid_bounds or self.bounds
        
        native_centre_coords = [(native_grid_bounds[dim_index] + native_grid_bounds[dim_index+2]) / 2.0 for dim_index in range(2)]
        utm_wkt = get_utm_wkt(native_centre_coords, self.wkt)
        
        return self.grid_points(grid_resolution=utm_grid_resolution, 
                                variables=variables,
                                native_grid_bounds=native_grid_bounds, 
                                resampling_method=resampling_method, 
                                grid_wkt=utm_wkt, 
                                point_step=point_step
                                )


    def utm_coords(self, coordinate_array, wkt=None):
        '''
        Function to convert coordinates to the appropriate UTM CRS
        @param coordinate_array: Array of shape (n, 2) or iterable containing coordinate pairs
        @param wkt: WKT for source CRS - default to native
       
        @return wkt: WKT for UTM CRS - default to native
        @return coordinate_array: Array of shape (n, 2) containing UTM coordinate pairs 
        '''
        wkt = wkt or self.wkt
        return utm_coords(coordinate_array, wkt)
    
    
    def coords2metres(self, coordinate_array, wkt=None):
        '''
        Function to calculate cumulative distance in metres from coordinates in specified CRS
        @param coordinate_array: Array of shape (n, 2) or iterable containing coordinate pairs
        @param wkt: WKT for coordinate CRS - default to native
        
        @return distance_array: Array of shape (n) containing cumulative distances from first coord
        '''
        wkt = wkt or self.wkt # Default to native CRS for coordinates

        _utm_wkt, utm_coord_array = utm_coords(coordinate_array, wkt)
        return coords2distance(utm_coord_array)


    def nearest_neighbours(self, coordinates, 
                           wkt=None, 
                           points_required=1, 
                           max_distance=None, 
                           secondary_mask=None):
        '''
        Function to determine nearest neighbours using cKDTree
        N.B: All distances are expressed in the native dataset CRS
        
        @param coordinates: two-element XY coordinate tuple, list or array
        @param wkt: Well-known text of coordinate CRS - defaults to native dataset CRS
        @param points_required: Number of points to retrieve. Default=1
        @param max_distance: Maximum distance to search from target coordinate - 
            STRONGLY ADVISED TO SPECIFY SENSIBLE VALUE OF max_distance TO LIMIT SEARCH AREA
        @param secondary_mask: Boolean array of same shape as point array used to filter points. None = no filter.
        
        @return distances: distances from the target coordinate for each of the points_required nearest points
        @return indices: point indices for each of the points_required nearest points
        '''
        if wkt:
            reprojected_coords = transform_coords(coordinates, wkt, self.wkt)
        else:
            reprojected_coords = coordinates
            
        if secondary_mask is None:
            secondary_mask = np.ones(shape=(self.point_count,), dtype=bool)
        else:
            assert secondary_mask.shape == (self.point_count,)        

        if max_distance:
            logger.debug('Computing spatial subset mask...')
            spatial_mask = self.get_spatial_mask([reprojected_coords[0] - max_distance,
                                                  reprojected_coords[1] - max_distance,
                                                  reprojected_coords[0] + max_distance,
                                                  reprojected_coords[1] + max_distance
                                                  ]
                                                 )
            
            point_indices = np.where(np.logical_and(spatial_mask,
                                                    secondary_mask
                                                    )
                                     )[0]
                                     
            if not len(point_indices):
                logger.debug('No points within distance {} of {}'.format(max_distance, reprojected_coords))
                return [], []
            
            # Set up KDTree for nearest neighbour queries
            logger.debug('Indexing spatial subset with {} points into KDTree...'.format(np.count_nonzero(spatial_mask)))
            kdtree = cKDTree(data=self.xycoords[point_indices])
            logger.debug('Finished indexing spatial subset into KDTree.')
        else:
            max_distance = np.inf
            if not self.kdtree:
                logger.debug('Indexing full dataset with {} points into KDTree...'.format(self.xycoords.shape[0]))
                self.kdtree = cKDTree(data=self.xycoords[np.where(secondary_mask)[0]])
                logger.debug('Finished indexing full dataset into KDTree.')
            kdtree = self.kdtree

            
        distances, indices = kdtree.query(x=np.array(reprojected_coords),
                                          k=points_required,
                                          distance_upper_bound=max_distance)
        
        if max_distance == np.inf:
            return distances, indices
        else: # Return indices of complete coordinate array, not the spatial subset
            return distances, np.where(spatial_mask)[0][indices]
            

    def get_lookup_mask(self, 
                        lookup_value_list, 
                        lookup_variable_name='line',
                        indexing_variable_name=None,
                        indexing_dimension='point'
                        ): 
        '''
        Function to return mask array based on lookup variable
        '''
        if lookup_variable_name:
            lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
            
            if lookup_variable.shape == (): # Scalar
                dimension = self.netcdf_dataset.get(indexing_dimension)
                assert dimension, 'Invalid indexing_dimension {} specified'.format(indexing_dimension)
                # Repeat boolean value across dimension size
                return np.array([lookup_variable[:] in lookup_value_list] * dimension.size)
        
            indexing_variable_name = indexing_variable_name or lookup_variable_name + '_index'
            
            try:
                indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
            except:
                raise BaseException('indexing_variable_name not supplied and cannot be inferred')
            
        elif indexing_variable_name:
            indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
            
            if hasattr(indexing_variable, 'lookup'): 
                # Get lookup variable name from variable attribute
                lookup_variable_name = indexing_variable.lookup
            elif indexing_variable_name.endswith('_index'):
                # Infer lookup variable name from indexing variable name
                lookup_variable_name = re.sub('_index$', '', indexing_variable_name)
            else:
                raise BaseException('lookup_variable_name not supplied and cannot be inferred')
            
            lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
        else:
            raise BaseException('Must supply either lookup_variable_name or indexing_variable_name')
        
        lookup_indices = np.arange(lookup_variable.shape[0])[np.in1d(lookup_variable[:], np.array(lookup_value_list))]
        logger.debug('lookup_indices: {}'.format(lookup_indices))  
          
        lookup_mask = np.in1d(indexing_variable, lookup_indices) 
        logger.debug('lookup_mask: {}'.format(lookup_mask))  
        return lookup_mask
                       

#===============================================================================
#     def lookup_mask_generator(self, 
#                         lookup_value_list, 
#                         lookup_variable_name='line',
#                         indexing_variable_name=None
#                         ): 
#         '''
#         Generator to yield mask array based on lookup variable for each of a list of lookup values
#         '''
#         if lookup_variable_name:
#             indexing_variable_name = indexing_variable_name or lookup_variable_name + '_index'
#             
#             try:
#                 indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
#             except:
#                 raise BaseException('indexing_variable_name not supplied and cannot be inferred')
#             
#         elif indexing_variable_name:
#             indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
#             
#             if hasattr(indexing_variable, 'lookup'): 
#                 # Get lookup variable name from variable attribute
#                 lookup_variable_name = indexing_variable.lookup
#             elif indexing_variable_name.endswith('_index'):
#                 # Infer lookup variable name from indexing variable name
#                 lookup_variable_name = re.sub('_index$', '', indexing_variable_name)
#             else:
#                 raise BaseException('lookup_variable_name not supplied and cannot be inferred')
#             
#         else:
#             raise BaseException('Must supply either lookup_variable_name or indexing_variable_name')
# 
#         lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
#         
#         for lookup_value in lookup_value_list:
#             lookup_indices = np.where(lookup_variable[:] == lookup_value)[0]
#             logger.debug('lookup_indices: {}'.format(lookup_indices))  
#               
#             lookup_mask = np.in1d(indexing_variable, lookup_indices) 
#             logger.debug('lookup_mask: {}'.format(lookup_mask))  
#             yield lookup_mask
#                        
#===============================================================================
    def get_index_mask(self, 
                        lookup_value_list, 
                        lookup_variable_name='line',
                        start_index_variable_name=None,
                        count_variable_name=None,
                        point_count=None
                        ): 
        '''
        Function to return mask array based on index variable
        '''
        try:
            lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
        except:
            raise BaseException('Invalid lookup_variable_name')

        start_index_variable_name = start_index_variable_name or lookup_variable_name + '_start_index'            
        try:
            start_index_variable = self.netcdf_dataset.variables[start_index_variable_name]
        except:
            raise BaseException('start_index_variable_name not supplied and cannot be inferred')
        
        count_variable_name = count_variable_name or lookup_variable_name + '_count'            
        try:
            count_variable = self.netcdf_dataset.variables[count_variable_name]
        except:
            raise BaseException('count_variable_name not supplied and cannot be inferred')

        point_count = point_count or self.netcdf_dataset.dimension['point'].size         

        
        lookup_indices = np.arange(lookup_variable.shape[0])[np.in1d(lookup_variable[:], lookup_value_list)]
        logger.debug('lookup_indices: {}'.format(lookup_indices))  
        start_indices = start_index_variable[lookup_indices]
        logger.debug('start_indices: {}'.format(start_indices))  
        counts = count_variable[lookup_indices]
        logger.debug('counts: {}'.format(counts))  
        
        # Build mask
        index_mask = np.zeros(shape=(point_count,), dtype='bool')  
        for lookup_index in range(len(lookup_indices)):
            index_mask[start_indices[lookup_index]:start_indices[lookup_index]+counts[lookup_index]] = True
            
        return index_mask 
    

    def expand_indexing_variable(self, 
                                 indexing_variable_name, 
                                 start_index=0, 
                                 end_index=0, 
                                 indexing_dimension='point',
                                 mask=None):
        '''
        Helper function to expand indexing variables and return an array of the required size
        '''       
        end_index = end_index or self.nc_dataset.dimensions[indexing_dimension].size
        index_range = end_index - start_index
         
        if mask is None: # No mask defined - take all points in range
            subset_mask = np.ones(shape=(index_range,), dtype='bool')
        else:
            subset_mask = mask[start_index:end_index]
        
        value_variable = self.nc_dataset.variables[indexing_variable_name]
        start_variable = self.nc_dataset.variables['{}_start_index'.format(indexing_variable_name)]
        count_variable = self.nc_dataset.variables['{}_point_count'.format(indexing_variable_name)]
            
        expanded_array = np.zeros(shape=(index_range,), 
                                  dtype=value_variable.dtype)
         
        # Assume monotonically increasing start indices to find all relevant indices
        indices = np.where(np.logical_and((start_variable[:] >= start_index),
                                          (start_variable[:] <= end_index)))[0]
         
        #logger.debug('indices: {}'.format(indices))
        for index in indices:
            start_index = max(start_variable[index]-start_index, 0)
            end_index = min(start_index+count_variable[index], index_range)
            expanded_array[start_index:end_index] = value_variable[index]
             
        #logger.debug('expanded_array: {}'.format(expanded_array))
        return expanded_array[subset_mask]
     
     
    def expand_lookup_variable(self, 
                               lookup_variable_name='line',
                               indexing_variable_name=None, 
                               start_index=0, 
                               end_index=0, 
                               mask=None,
                               indexing_dimension='point'):
        '''
        Function to expand lookup variables and return an array of the required size
        '''
        if lookup_variable_name:
            lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
            
            if lookup_variable.shape == (): # Scalar
                dimension = self.netcdf_dataset.get(indexing_dimension)
                assert dimension, 'Invalid indexing_dimension {} specified'.format(indexing_dimension)
                # Repeat boolean value across dimension size
                return np.array([lookup_variable[:]] * dimension.size)
        
            indexing_variable_name = indexing_variable_name or lookup_variable_name + '_index'
            
            try:
                indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
            except:
                raise BaseException('indexing_variable_name not supplied and cannot be inferred')
            
        elif indexing_variable_name:
            indexing_variable = self.netcdf_dataset.variables[indexing_variable_name]
            
            if hasattr(indexing_variable, 'lookup'): 
                # Get lookup variable name from variable attribute
                lookup_variable_name = indexing_variable.lookup
            elif indexing_variable_name.endswith('_index'):
                # Infer lookup variable name from indexing variable name
                lookup_variable_name = re.sub('_index$', '', indexing_variable_name)
            else:
                raise BaseException('lookup_variable_name not supplied and cannot be inferred')
            
            lookup_variable = self.netcdf_dataset.variables[lookup_variable_name]
        else:
            raise BaseException('Must supply either lookup_variable_name or indexing_variable_name')
             
        end_index = end_index or indexing_variable.shape[0] # Usually this will be the point count
        index_range = end_index - start_index
         
        if mask is None: # No mask defined - take all points in range
            subset_mask = np.ones(shape=(index_range,), dtype='bool')
        else:
            subset_mask = mask[start_index:end_index]
        
        return lookup_variable[indexing_variable[start_index:end_index][subset_mask]]

                       
#===============================================================================
#     def read_points(self, start_index, end_index, mask=None):
#         '''
#         '''
#     
#     
#         # Start of read_points function
#         self.index_range = end_index - start_index
#         
#         if mask is None: # No mask defined - take all points in range
#             subset_mask = np.ones(shape=(self.index_range,), dtype='bool')
#         else:
#             subset_mask = mask[start_index:end_index]
#             self.index_range = np.count_nonzero(subset_mask)
#             
#         # If no points to retrieve, don't read anything
#         if not self.index_range:
#             logger.debug('No points to retrieve - all masked out')            
#             return
# 
#         for field_definition in self.field_definitions:
#             short_name = field_definition['short_name']
#             variable_name = field_definition['variable_name']
#             #logger.debug('variable_name: {}'.format(variable_name))
#             
#             variable = self.nc_dataset.variables[variable_name]
#             
#             if len(variable.shape) == 0: # Scalar variable - repeat value for each point
#                 data = variable[:]
#                 self.cache[short_name] = np.array([data] * self.index_range)
#                 
#             elif len(variable.shape) in [1, 2]: # 1D or 2D variable
#                 # Indexing array variable
#                 if ('point' not in variable.dimensions) and (short_name in self.settings['index_fields']): 
#                     self.cache[short_name] = expand_indexing_variable(variable_name)[subset_mask]
#                     
#                 # Lookup array variable
#                 elif ('point' in variable.dimensions) and (short_name in self.settings['lookup_fields'] or variable_name.endswith('_index')): 
#                     self.cache[short_name] = expand_lookup_variable(variable_name)[subset_mask]
# 
#                 # A data array variable
#                 elif ('point' in variable.dimensions) and (variable_name not in self.settings['index_fields']):
#                     data = variable[start_index:end_index]
#                     
#                     # Include fill_values if array is masked
#                     if type(data) == np.ma.core.MaskedArray:
#                         data = data.data
#                         
#                     self.cache[short_name] = data[subset_mask]
#             else:
#                 raise BaseException('Invalid dimensionality for variable {}'.format(variable_name))     
#         
#         #logger.debug('self.cache: {}'.format(pformat(self.cache)))
# 
#     def point_generator(self, mask=None):
#         '''
#         '''
#         pass
#===============================================================================
        
    