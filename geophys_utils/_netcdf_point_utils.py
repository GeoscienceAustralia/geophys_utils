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
import numexpr as ne
import math
import os
import sys
import re
import tempfile
from collections import OrderedDict
from pprint import pformat
from scipy.interpolate import griddata
from geophys_utils._crs_utils import transform_coords, get_utm_wkt, get_reprojected_bounds, get_spatial_ref_from_wkt
from geophys_utils._transect_utils import utm_coords, coords2distance
from geophys_utils._netcdf_utils import NetCDFUtils, METADATA_CRS
from geophys_utils._polygon_utils import points2convex_hull
from geophys_utils._concave_hull import concaveHull
from shapely.geometry import shape
from scipy.spatial.ckdtree import cKDTree
from shapely.geometry import Polygon, MultiPoint, MultiPolygon
from shapely.geometry.polygon import asPolygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform
import shapely.wkt
import logging

# Setup logging handlers if required
logger = logging.getLogger(__name__) # Get logger
logger.setLevel(logging.INFO) # Initial logging level for this module

try:
    import memcache
except ImportError:
    logger.debug('Unable to import memcache. AWS-specific functionality will not be enabled')
    memcache = None


# Default number of points to read per chunk when retrieving data
DEFAULT_READ_CHUNK_SIZE = 8192

# Set this to a number other than zero for testing
POINT_LIMIT = 0
    
# Metadata shape generation parameters
SHAPE_BUFFER_DISTANCE = 0.02 # Distance to buffer (kerf) shape out then in again (in degrees)
SHAPE_OFFSET = 0.0005 # Distance to buffer (kerf) final shape outwards (in degrees)
SHAPE_SIMPLIFY_TOLERANCE = 0.0005 # Length of shortest line in shape (in degrees)
SHAPE_MAX_POLYGONS=5
SHAPE_MAX_VERTICES=1000
SHAPE_ORDINATE_DECIMAL_PLACES = 6 # Number of decimal places for shape vertex ordinates


class NetCDFPointUtils(NetCDFUtils):
    '''
    NetCDFPointUtils class to do various fiddly things with NetCDF geophysics point data files.
    '''
    CACHE_VARIABLE_PARAMETERS = {'complevel': 4,
                                   'zlib': True,
                                   'fletcher32': True,
                                   'shuffle': True,
                                   'endian': 'little',
                                   }

    def __init__(self, 
                 netcdf_dataset,
                 memcached_connection=None,
                 enable_disk_cache=None,
                 enable_memory_cache=True,
                 cache_path=None,
                 s3_bucket=None,
                 debug=False):
        '''
        NetCDFPointUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a point dataset
        @parameter enable_disk_cache: Boolean parameter indicating whether local cache file should be used, or None for default 
        @parameter enable_memory_cache: Boolean parameter indicating whether values should be cached in memory or not.
        @parameter debug: Boolean parameter indicating whether debug output should be turned on or not
        '''
        # Start of init function - Call inherited constructor first
        super().__init__(netcdf_dataset=netcdf_dataset, 
                         debug=debug
                         )
        
        logger.debug('Running NetCDFPointUtils constructor')
        
        if memcache is not None:
            self.memcached_connection = memcached_connection
        else:
            self.memcached_connection = None

        self.s3_bucket = s3_bucket

        self.cache_path = cache_path or os.path.join(os.path.join(tempfile.gettempdir(), 'NetCDFPointUtils'),
                                                     re.sub('\W', '_', os.path.splitext(self.nc_path)[0])) + '_cache.nc'

        self.cache_basename = os.path.join(self.cache_path,
                                           re.sub('\W', '_', os.path.splitext(self.nc_path)[0]))


        #logger.debug('self.cache_path')
        #logger.debug(self.cache_path)
        #logger.debug('self.cache_path: {}'.format(self.cache_path))

        self.enable_memory_cache = enable_memory_cache
        
        # If caching is not explicitly specified, enable it for OPeNDAP access
        if enable_disk_cache is None:
            self.enable_disk_cache = self.opendap
        else:
            self.enable_disk_cache = enable_disk_cache

        # Initialise private property variables to None until set by property getter methods
        self._xycoords = None         
        self._point_variables = None        
        self._data_variable_list = None
        self._kdtree = None

        # Determine exact spatial bounds
        xycoords = self.xycoords
        xmin = np.nanmin(xycoords[:,0])
        xmax = np.nanmax(xycoords[:,0])
        ymin = np.nanmin(xycoords[:,1])
        ymax = np.nanmax(xycoords[:,1])

        # Create nested list of bounding box corner coordinates
        self.native_bbox = [[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]]

        # Define bounds
        self.bounds = [xmin, ymin, xmax, ymax]
            
        self.point_count = self.netcdf_dataset.dimensions['point'].size
               
        
    #===========================================================================
    # def __del__(self):
    #     '''
    #     NetCDFPointUtils Destructor
    #     '''
    #     if self.enable_disk_cache:
    #         try:
    #             cache_file_path = self._nc_cache_dataset.filepath()
    #             self._nc_cache_dataset.close()
    #             os.remove(cache_file_path)
    #         except:
    #             pass
    #===========================================================================
        
    def fetch_array(self, source_variable, dest_array=None):
        '''
        Helper function to retrieve entire 1D array in pieces < self.max_bytes in size
        @param source_variable: netCDF variable from which to retrieve data
        '''
        source_len = source_variable.shape[0]
        pieces_required = int(math.ceil((source_variable[0].itemsize * source_len) / self.max_bytes))
        max_elements = source_len // pieces_required
        
        # Reduce max_elements to fit within chunk boundaries if possible
        if pieces_required > 1 and hasattr(source_variable, '_ChunkSizes'):
            chunk_size = (source_variable._ChunkSizes 
                          if type(source_variable._ChunkSizes) in [int, np.int32] 
                          else source_variable._ChunkSizes[0]
                          )
            chunk_count = max(max_elements // chunk_size, 
                              1)
            max_elements = min(chunk_count * chunk_size, 
                               max_elements)
            pieces_required = int(math.ceil(source_len / max_elements))
        
        logger.debug('Fetching {} pieces containing up to {} {} array elements.'.format(pieces_required, max_elements, source_variable.name))
        
        if dest_array is None:
            dest_array = np.zeros((source_len,), dtype=source_variable.dtype)

        # Copy array in pieces
        start_index = 0
        while start_index < source_len:
            end_index = min(start_index + max_elements, source_len)
            logger.debug('Retrieving {} array elements {}:{}'.format(source_variable.name, start_index, end_index))
            array_slice = slice(start_index, end_index)
            dest_array[array_slice] = source_variable[array_slice]
            start_index += max_elements
            
        return dest_array
        
    def get_polygon(self):
        '''
        Returns GML representation of convex hull polygon for dataset
        '''
        return 'POLYGON((' + ', '.join([' '.join(
            ['%.4f' % ordinate for ordinate in coordinates]) 
            for coordinates in self.get_convex_hull()]) + '))'
        
    
    def get_spatial_mask(self, bounds, bounds_wkt=None):
        '''
        Return boolean mask of dimension 'point' for all coordinates within specified bounds and CRS
        @parameter bounds: Either an iterable containing [<xmin>, <ymin>, <xmax>, <ymax>] or a shapely (multi)polygon
        @parameter bounds_wkt: WKT for bounds CRS. Defaults to dataset native CRS
        :return mask: Boolean array of size n
        '''
        
        #TODO: Deal with this in a more high-level way
        POINT_CHUNK_SIZE = 1048576 # Number of points to check at any one time to keep memory usage down
            
        def get_intersection_mask(points, geometry):
            """
            Determine if points lie inside (multi)polygon
            N.B: points and geometry must be in the same CRS
            :param points: 2 x n array of input coordinates
            :param geometry: (multi)polygon
            :return mask: Boolean array of size n 
            """
            mask = np.zeros(shape=(points.shape[0]), dtype=np.bool)
            
            chunk_start_index = 0
            while chunk_start_index < len(points):
                chunk_end_index = min(chunk_start_index + POINT_CHUNK_SIZE, len(points))
                logger.debug('Checking spatial containment for points {} to {} of {}'.format(chunk_start_index, chunk_end_index-1, len(points)))
                intersection_points = np.array(MultiPoint(points[slice(chunk_start_index, chunk_end_index)]).intersection(geometry))
                
                #TODO: Find out if there's a better way of getting the mask from the intersection points
                # Note that this method would have some issues with duplicated coordinates, but there shouldn't be any
                logger.debug('Computing partial mask from {} intersection points'.format(len(intersection_points)))
                _x_values, x_indices, _x_intersection_indices = np.intersect1d(points.flatten()[0::2], intersection_points.flatten()[0::2], return_indices=True)
                _y_values, y_indices, _y_intersection_indices = np.intersect1d(points.flatten()[1::2], intersection_points.flatten()[1::2], return_indices=True)
                intersection_indices = np.intersect1d(x_indices, y_indices, return_indices=False)
                mask[intersection_indices] = True
                    
                chunk_start_index = chunk_end_index

            return mask
            
            
        coordinates = self.xycoords # Don't transform these - do all spatial operations in native CRS
        #logger.debug('coordinates = {}'.format(coordinates))
    
    
        if isinstance(bounds, BaseGeometry): # Process shapely (multi)polygon bounds    
            if bounds_wkt is None:
                native_crs_bounds = bounds
            else:
                logger.debug('Original bounds = {}'.format(bounds))
                native_crs_bounds = transform((lambda x, y: transform_coords([x, y], bounds_wkt, self.wkt)), 
                                              bounds)
                
            logger.debug('native_crs_bounds = {}'.format(native_crs_bounds))

            # Shortcut the whole process if the extents are within the bounds geometry       
            if asPolygon(self.native_bbox).within(native_crs_bounds):
                logger.debug('Dataset is completely contained within bounds')
                return np.ones(shape=(len(coordinates),), dtype=np.bool)
                
            bounds_half_size = abs(np.array([native_crs_bounds.bounds[2] - native_crs_bounds.bounds[0], 
                                             native_crs_bounds.bounds[3] - native_crs_bounds.bounds[1]])) / 2.0
            bounds_centroid = np.array(native_crs_bounds.centroid.coords[0])
            #logger.debug('bounds_half_size = {}, bounds_centroid = {}'.format(bounds_half_size, bounds_centroid))
            
            # Limit the points checked to those within the same rectangular extent (for speed)
            # Set mask element to true for each point which is <= bounds_half_size distance from bounds_centroid
            mask = np.all(ne.evaluate("abs(coordinates - bounds_centroid) <= bounds_half_size"), axis=1)
            logger.debug('{}/{} points found in initial bounding box intersection'.format(np.count_nonzero(mask), len(coordinates)))
            
            # Apply sub-mask for all points within bounds geometry
            (mask[mask])[~get_intersection_mask(coordinates[mask], native_crs_bounds)] = False
            #logger.debug('Final shape mask = {}'.format(mask))
            
        else: # Process four-element bounds iterable if possible
            assert len(bounds) == 4, 'Invalid bounds iterable: {}. Must be of form [<xmin>, <ymin>, <xmax>, <ymax>]'.format(bounds)
            
            native_crs_bounds = transform_coords(np.array(bounds).reshape((2,2)), bounds_wkt, self.wkt).reshape((4, 1)) # Transform as [xmin, ymin], [xmax, ymax]]
                
            if (self.bounds[0] >= native_crs_bounds[0]
                and self.bounds[1] >= native_crs_bounds[1]
                and self.bounds[2] <= native_crs_bounds[2]
                and self.bounds[3] <= native_crs_bounds[3]
                ):
                logger.debug('Dataset is completely contained within bounds')
                return np.ones(shape=(len(coordinates),), dtype=np.bool)
                
        
            bounds_half_size = abs(np.array([native_crs_bounds[2] - native_crs_bounds[0], native_crs_bounds[3] - native_crs_bounds[1]])) / 2.0
            bounds_centroid = np.array([native_crs_bounds[0], native_crs_bounds[1]]) + bounds_half_size
            
            # Return true for each point which is <= bounds_half_size distance from bounds_centroid
            mask = np.all(ne.evaluate("abs(coordinates - bounds_centroid) <= bounds_half_size"), axis=1)  
                      
        logger.debug('{}/{} points found in final mask'.format(np.count_nonzero(mask), len(coordinates)))
        return mask
        
        
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
            reprojected_grid_bounds = get_reprojected_bounds(native_grid_bounds, self.wkt, grid_wkt)
        elif reprojected_grid_bounds:
            native_grid_bounds = get_reprojected_bounds(reprojected_grid_bounds, grid_wkt, self.wkt)
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

        spatial_subset_mask = self.get_spatial_mask(get_reprojected_bounds(expanded_grid_bounds, grid_wkt, self.wkt))
        
        # Create grids of Y and X values. Note YX ordering and inverted Y
        # Note GRID_RESOLUTION/2.0 fudge to avoid truncation due to rounding error
        grid_y, grid_x = np.mgrid[pixel_centre_bounds[3]:pixel_centre_bounds[1]-grid_resolution/2.0:-grid_resolution, 
                                 pixel_centre_bounds[0]:pixel_centre_bounds[2]+grid_resolution/2.0:grid_resolution]

        # Skip points to reduce memory requirements
        #TODO: Implement function which grids spatial subsets.
        point_subset_mask = np.zeros(shape=(self.netcdf_dataset.dimensions['point'].size,), dtype=bool)
        point_subset_mask[0:-1:point_step] = True
        point_subset_mask = np.logical_and(spatial_subset_mask, point_subset_mask)
        
        coordinates = self.xycoords[point_subset_mask]
        # Reproject coordinates if required
        if grid_wkt is not None:
            # N.B: Be careful about XY vs YX coordinate order         
            coordinates = np.array(transform_coords(coordinates[:], self.wkt, grid_wkt))

        # Interpolate required values to the grid - Note YX ordering for image
        grids = {}
        for variable in [self.netcdf_dataset.variables[var_name] for var_name in variables]:
            grids[variable.name] = griddata(coordinates[:,::-1],
                                  variable[:][point_subset_mask], #TODO: Check why this is faster than direct indexing
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


    def get_convex_hull(self, to_wkt=None):
        '''
        Function to return vertex coordinates of a convex hull polygon around all points
        Implements abstract base function in NetCDFUtils 
        @param to_wkt: CRS WKT for shape
        '''
        return points2convex_hull(transform_coords(self.xycoords, self.wkt, to_wkt))
    
    
    def get_concave_hull_old(self, to_wkt=None, smoothness=None):
        """\
        Returns the concave hull (as a shapely polygon) of all points.
        Implements abstract base function in NetCDFUtils
        @param to_wkt: CRS WKT for shape
        @param smoothness: distance to buffer (kerf) initial shape outwards then inwards to simplify it
        """
        hull = concaveHull(transform_coords(self.xycoords, self.wkt, to_wkt))
        result = shape({'type': 'Polygon', 'coordinates': [hull.tolist()]})

        if smoothness is None:
            return result

        return Polygon(result.buffer(smoothness).exterior).buffer(-smoothness)



    def get_concave_hull(self, to_wkt=None, buffer_distance=0.02, offset=0.0005, tolerance=0.0005, cap_style=1,
                         join_style=1, max_polygons=5, max_vertices=1000):
        """\
        Returns the concave hull (as a shapely polygon) of points with data.
        Implements abstract base function in NetCDFUtils
        @param to_wkt: CRS WKT for shape
        @param buffer_distance: distance to buffer (kerf) initial shape outwards then inwards to simplify it
        @param offset: Final offset of final shape from original lines
        @param tolerance: tolerance for simplification
        @param cap_style: cap_style for buffering. Defaults to round
        @param join_style: join_style for buffering. Defaults to round
        @param max_polygons: Maximum number of polygons to accept. Will keep doubling buffer_distance until under this limit. 0=unlimited.
        @param max_vertices: Maximum number of vertices to accept. Will keep doubling buffer_distance until under this limit. 0=unlimited.
        @return shapely.geometry.shape: Geometry of concave hull
        """
        assert not max_polygons or buffer_distance > 0, 'buffer_distance must be greater than zero if number of polygons is limited'  # Avoid endless recursion

        def get_offset_geometry(geometry, buffer_distance, offset, tolerance, cap_style, join_style, max_polygons,
                                max_vertices):
            '''\
            Helper function to return offset geometry. Will keep trying larger buffer_distance values until there is a manageable number of polygons
            '''
            logger.debug('Computing offset geometry with buffer_distance = {}'.format(buffer_distance))

            #transform_coords(line_vertices, self.wkt, to_wkt))
            offset_geometry = geometry.buffer(buffer_distance, cap_style=cap_style, join_style=join_style).simplify(
                tolerance)
            offset_geometry = offset_geometry.buffer(offset - buffer_distance, cap_style=cap_style,
                                                     join_style=join_style).simplify(tolerance)

            # Discard any internal polygons
            if type(offset_geometry) == MultiPolygon:
                polygon_list = []
                for polygon in offset_geometry:
                    polygon = Polygon(polygon.exterior)
                    polygon_is_contained = False
                    for list_polygon in polygon_list:
                        polygon_is_contained = list_polygon.contains(polygon)
                        if polygon_is_contained:
                            break
                        elif polygon.contains(list_polygon):
                            polygon_list.remove(list_polygon)
                            break

                    if not polygon_is_contained:
                        polygon_list.append(polygon)

                if len(polygon_list) == 1:
                    offset_geometry = polygon_list[0]  # Single polygon
                else:
                    offset_geometry = MultiPolygon(polygon_list)

            elif type(offset_geometry) == Polygon:
                offset_geometry = Polygon(offset_geometry.exterior)
            else:
                raise ValueError('Unexpected type of geometry: {}'.format(type(offset_geometry)))

            # Keep doubling the buffer distance if there are too many polygons
            if (
                    (max_polygons and type(offset_geometry) == MultiPolygon and len(offset_geometry) > max_polygons)
                    or
                    (max_vertices and type(offset_geometry) == MultiPolygon and
                     sum([len(polygon.exterior.coords)
                          # + sum([len(interior_ring.coords) for interior_ring in polygon.interiors])
                          for polygon in offset_geometry]) > max_vertices)
                    or
                    (max_vertices and type(offset_geometry) == Polygon and
                     (len(offset_geometry.exterior.coords)
                     # + sum([len(interior_ring.coords) for interior_ring in offset_geometry.interiors])
                     ) > max_vertices)
            ):
                return get_offset_geometry(geometry, buffer_distance * 2, offset, tolerance, cap_style, join_style,
                                           max_polygons, max_vertices)

            return offset_geometry

        return get_offset_geometry(Polygon(self._xycoords), buffer_distance, offset, tolerance, cap_style, join_style,
                                   max_polygons, max_vertices)

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

        if max_distance: # max_distance has been specified
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
        else: # Consider ALL points
            max_distance = np.inf
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
            
            if (lookup_variable.shape == () 
                or ((len(lookup_variable.shape) == 1) and (lookup_variable.dtype == '|S1'))): # Scalar or string array
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
        
        # Handle special case for string arrays via OPeNDAP
        if self.opendap and (lookup_variable.dtype == 'S1') and (len(lookup_variable.shape) == 2):
            # Convert 2D byte array into 1D array of unicode strings - needed for OPeNDAP
            lookup_array = np.array([bytestring[bytestring != b''].tostring().decode('UTF8') for bytestring in lookup_variable[:]])
            # OPeNDAP will truncate strings to 64 characters - truncate search strings to match
            lookup_indices = np.arange(lookup_array.shape[0])[np.in1d(lookup_array, np.array([lookup_value[0:64] 
                                                                                              for lookup_value in lookup_value_list]))]
        else:
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

        point_count = point_count or self.netcdf_dataset.dimensions['point'].size         

        
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
                dimension = self.netcdf_dataset.dimensions.get(indexing_dimension)
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
            
        lookup_variable.set_auto_mask(False)
        indexing_variable.set_auto_mask(False)
            
        result_array = lookup_variable[:][indexing_variable[start_index:end_index][subset_mask]] # Need to index numpy array, not netCDF variable

        # Convert 2D byte array into 1D array of unicode strings - needed for OPeNDAP
        if result_array.dtype == 'S1':
            result_array = np.array([bytestring[bytestring != b''].tostring().decode('UTF8') for bytestring in result_array])
            
        return result_array
                       
    def chunk_point_data_generator(self, 
                                   start_index=0, 
                                   end_index=0,
                                   field_list=None,
                                   mask=None,
                                   yield_variable_attributes_first=False):
        '''
        Generator to optionally yield variable attributes followed by all point data for the specified point index range
        Used to retrieve data as chunks for outputting as point-wise lists of lists
        @param start_index: start point index of range to read 
        @param end_index: end point index of range to read. Defaults to number of points
        @param field_list: Optional list of field names to read. Default is None for all variables 
        @param mask: Optional Boolean mask array to subset points
        @param yield_variable_attributes_first: Boolean flag to determine whether variable attribute dict is yielded first. Defaults to False
        
        @yield variable_attributes: dict of netCDF variable attributes. Optionally the first item yielded if yield_variable_attributes_first is True
        @yield point_value_list: List of single values for 1D variables or sub-lists for 2D variables for a single point
        '''
        # Start of point_data_generator function
        end_index = end_index or self.point_count
        index_range = end_index - start_index
         
        if mask is None: # No mask defined - take all points in range
            subset_mask = np.ones(shape=(index_range,), dtype='bool')
        else:
            subset_mask = mask[start_index:end_index]
            index_range = np.count_nonzero(subset_mask)
             
        # If no points to retrieve, don't read anything
        if not index_range:
            logger.debug('No points to retrieve for point indices {}-{}: All masked out'.format(start_index, end_index-1))            
            return
        
        # Generate full field list if None provided
        if not field_list:
            field_list = [variable.name 
                          for variable in self.netcdf_dataset.variables.values()
                          if (not len(variable.dimensions) # Scalar variable
                              or variable.dimensions[0] == 'point' # Variable is of point dimension
                              or (variable.dimensions[0] + '_index' in self.netcdf_dataset.variables.keys() # Variable has an index variable
                                  and len(self.netcdf_dataset.variables[variable.dimensions[0] + '_index'].dimensions) # index variable is not a scalar
                                  and self.netcdf_dataset.variables[variable.dimensions[0] + '_index'].dimensions[0] == 'point' # index variable is of point dimension
                                  )
                              )
                          and not variable.name.endswith('_index') 
                            and not hasattr(variable, 'lookup') # Variable is not an index variable
                          and not variable.name in NetCDFUtils.CRS_VARIABLE_NAMES 
                            and not re.match('ga_.+_metadata', variable.name) # Not an excluded variable
                          ]
 
        logger.debug('field_list: {}'.format(field_list))
        
        variable_attributes = OrderedDict()
        memory_cache = OrderedDict()
        for variable_name in field_list:
            variable = self.netcdf_dataset.variables.get(variable_name)
            if variable is None:
                logger.warning('Variable {} does not exist. Skipping.'.format(variable_name))
                continue
            #logger.debug('variable_name: {}'.format(variable_name))
            
            # Scalar variable
            if len(variable.shape) == 0:
                # Skip CRS variable
                if variable_name in NetCDFUtils.CRS_VARIABLE_NAMES or re.match('ga_.+_metadata', variable_name):
                    continue 
                
                # Repeat scalar value for each point
                data_array = variable[:]
                memory_cache[variable_name] = np.array([data_array] * index_range)
                 
            else: # nD array variable
                if (variable.dimensions[0] != 'point'): # Variable is NOT of point dimension - must be lookup
                    memory_cache[variable_name] = self.expand_lookup_variable(lookup_variable_name=variable_name, 
                                                                              start_index=start_index, 
                                                                              end_index=end_index, 
                                                                              mask=mask)                     
                else: # 'point' is in variable.dimensions - "normal" variable
                    data_array = variable[start_index:end_index]
                     
                    # Include fill_values if array is masked
                    if type(data_array) == np.ma.core.MaskedArray:
                        data_array = data_array.data
                         
                    memory_cache[variable_name] = data_array[subset_mask]
               
            if yield_variable_attributes_first:
                variable_attributes[variable_name] = dict(variable.__dict__)
            
        logger.debug('variable_attributes: {}'.format(pformat(variable_attributes)))
        logger.debug('memory_cache: {}'.format(pformat(memory_cache)))
        
        if yield_variable_attributes_first:
            yield variable_attributes
        
        for index in range(index_range):
            point_value_list = []
            for variable_name, variable in iter(memory_cache.items()):
                data_array = variable[index]
                
                # Convert array to string if required
                if type(data_array) == np.ndarray and data_array.dtype == object:
                    data_array = str(data_array)

                point_value_list.append(data_array)
                            
            yield point_value_list
            
        logger.debug('{} points read for point indices {}-{}'.format(index_range, start_index, end_index-1))        
        
 
    def all_point_data_generator(self,
                                 field_list=None,
                                 mask=None,
                                 read_chunk_size=None,
                                 yield_variable_attributes_first=True):
        '''
        Generator to yield variable attributes followed by lists of values for all points
        @param field_list: Optional list of field names to read. Default is None for all variables 
        @param mask: Optional Boolean mask array to subset points
        @param read_chunk_size: Number of points to read from the netCDF per chunk (for greater efficiency than single point reads)
        @param yield_variable_attributes_first: Boolean flag to determine whether variable attribute dict is yielded first. Defaults to True
        
        @yield variable_attributes: dict of netCDF variable attributes. Optionally the first item yielded if yield_variable_attributes_first is True
        @yield point_value_list: List of single values for 1D variables or sub-lists for 2D variables for a single point
        '''
        read_chunk_size = read_chunk_size or DEFAULT_READ_CHUNK_SIZE
        
        # Process all chunks
        point_count = 0
        for chunk_index in range(self.point_count // read_chunk_size + 1):
            for line in self.chunk_point_data_generator(field_list=field_list,
                                                        start_index=chunk_index*read_chunk_size,
                                                        end_index=min((chunk_index+1)*read_chunk_size,
                                                                      self.point_count
                                                                      ),
                                                        mask=mask,
                                                        yield_variable_attributes_first=yield_variable_attributes_first
                                             ):
                if not yield_variable_attributes_first:
                    point_count += 1
                
                yield_variable_attributes_first = False # Only yield variable attributes from the first chunk

                #logger.debug('line: {}'.format(line))
                yield line
            
                if POINT_LIMIT and (point_count >= POINT_LIMIT):
                    break                    
                
            if POINT_LIMIT and (point_count >= POINT_LIMIT):
                break
        
        logger.debug('{} points read from netCDF file {}'.format(point_count, self.nc_path))

    def get_xy_coord_values(self):
        '''
        Function to return a full in-memory coordinate array from source dataset
        '''
        logger.debug('Reading xy coordinates from source dataset')
            
        xycoord_values = np.zeros(shape=(len(self.x_variable), 2), dtype=self.x_variable.dtype)
        self.fetch_array(self.x_variable, xycoord_values[:,0])
        self.fetch_array(self.y_variable, xycoord_values[:,1])
        
        # Deal with netCDF4 Datasets that have had set_auto_mask(False) called
        if hasattr(self.x_variable, '_FillValue'):
            xycoord_values[:,0][xycoord_values[:,0] == self.x_variable._FillValue] = np.nan
        if hasattr(self.y_variable, '_FillValue'):
            xycoord_values[:,1][xycoord_values[:,1] == self.y_variable._FillValue] = np.nan
        
        return xycoord_values    

    @property
    def xycoords(self):
        '''
        Property getter function to return pointwise array of XY coordinates
        The order of priority for retrieval is memory, memcached, disk cache then dataset.
        '''
        xycoords = None
            # assert np.allclose(arr, arr_down)

        if self.enable_memory_cache and self._xycoords is not None:
            #logger.debug('Returning memory cached coordinates')
            return self._xycoords

        elif self.memcached_connection is not None:
            coord_cache_key = self.cache_basename + '_xycoords'
            
            logger.debug("hit xycoords propery code")
            logger.debug(self.memcached_connection)

            xycoords = self.memcached_connection.get(coord_cache_key)
            if xycoords is not None:
                # self.memcached_connection.get(self.cache_path) is True:
                logger.debug('memcached key found at {}'.format(coord_cache_key))
                #logger.debug('xycoords: {}'.format(xycoords))
            else:
                xycoords = self.get_xy_coord_values()
                logger.debug("key not found at {}. adding key and value".format(coord_cache_key))
                self.memcached_connection.add(coord_cache_key, xycoords)


        elif self.enable_disk_cache:
            if os.path.isfile(self.cache_path):
                # Cached coordinate file exists - read it
                cache_dataset = netCDF4.Dataset(self.cache_path, 'r')

                #assert cache_dataset.source == self.nc_path, 'Source mismatch: cache {} vs. dataset {}'.format(cache_dataset.source, self.nc_path)

                if 'xycoords' in cache_dataset.variables.keys():
                    xycoords = cache_dataset.variables['xycoords'][:]
                    logger.debug('Read {} coordinates from cache file {}'.format(xycoords.shape[0], self.cache_path))
                else:
                    logger.debug('Unable to read xycoords variable from netCDF cache file {}'.format(self.cache_path))
                cache_dataset.close()
            else:
                logger.debug('NetCDF cache file {} does not exist'.format(self.cache_path))

            if xycoords is None:
                xycoords = self.get_xy_coord_values() # read coords from source file

                os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
                if os.path.isfile(self.cache_path):
                    cache_dataset = netCDF4.Dataset(self.cache_path, 'r+')
                else:
                    cache_dataset = netCDF4.Dataset(self.cache_path, 'w')

                if not hasattr(cache_dataset, 'source'):
                    cache_dataset.source = self.nc_path

                #assert cache_dataset.source == self.nc_path, 'Source mismatch: cache {} vs. dataset {}'.format(cache_dataset.source, self.nc_path)

                if 'point' not in cache_dataset.dimensions.keys():
                    cache_dataset.createDimension(dimname='point', size=xycoords.shape[0])

                if 'xy' not in cache_dataset.dimensions.keys():
                    cache_dataset.createDimension(dimname='xy', size=xycoords.shape[1])

                if 'xycoords' not in cache_dataset.variables.keys():
                    cache_dataset.createVariable('xycoords',
                                                 xycoords.dtype,
                                                 dimensions=['point', 'xy'],
                                                 **self.CACHE_VARIABLE_PARAMETERS
                                                 )
                cache_dataset.variables['xycoords'][:] = xycoords # Write coords to cache file
                cache_dataset.close()
                logger.debug('Saved {} coordinates to cache file {}'.format(xycoords.shape[0], self.cache_path))
            
        else: # No caching - read coords from source file
            xycoords = self.get_xy_coord_values()

        if self.enable_memory_cache:
            self._xycoords = xycoords
            
        return xycoords
        
    @property
    def point_variables(self):
        '''
        Property getter function to return point_variables as required
        '''
        if not self._point_variables:
            logger.debug('Setting point_variables property')
            self._point_variables = list([var_name for var_name in self.netcdf_dataset.variables.keys() 
                                          if 'point' in self.netcdf_dataset.variables[var_name].dimensions
                                          and var_name not in ['latitude', 'longitude', 'easting', 'northing', 'point', 'fiducial', 'flag_linetype']
                                          ])
        return self._point_variables
         
         
    @property
    def data_variable_list(self):
        '''
        Property getter function to return data_variable_list as required
        '''
        if not self._data_variable_list:
            logger.debug('Setting data_variable_list property')
            self._data_variable_list = [key for key, value in self.netcdf_dataset.variables.items()
                                       if 'point' in value.dimensions]
        return self._data_variable_list

        
    @property
    def kdtree(self):
        '''
        Property getter function to return data_variable_list as required
        '''
        if not self._kdtree:
            logger.debug('Indexing full dataset with {} points into KDTree...'.format(self.xycoords.shape[0]))
            self._kdtree = cKDTree(data=self.xycoords, balanced_tree=False)
            logger.debug('Finished indexing full dataset into KDTree.')
        return self._kdtree

    def copy(self, 
             nc_out_path, 
             datatype_map_dict={},
             variable_options_dict={},
             dim_range_dict={},
             dim_mask_dict={},
             nc_format=None,
             limit_dim_size=False,
             var_list=[],
             empty_var_list=[],
             to_crs=None       
        ):
        '''
        Function to copy a netCDF dataset to another one with potential changes to size, format, 
            variable creation options and datatypes.
            
            @param nc_out_path: path to netCDF output file 
            @param to_crs: WKT of destination CRS

        '''  
        
        if var_list:
            expanded_var_list = list(set(
                var_list + 
                NetCDFUtils.X_DIM_VARIABLE_NAMES + 
                NetCDFUtils.Y_DIM_VARIABLE_NAMES +
                NetCDFUtils.CRS_VARIABLE_NAMES +
                ['line', 'line_index'] # Always include line numbers (This really should be in an overridden function in NetCDFLineUtils)
                ))
        else:
            expanded_var_list = var_list
        
        # Call inherited NetCDFUtils method
        super().copy( 
             nc_out_path, 
             datatype_map_dict=datatype_map_dict,
             variable_options_dict=variable_options_dict,
             dim_range_dict=dim_range_dict,
             dim_mask_dict=dim_mask_dict,
             nc_format=nc_format,
             limit_dim_size=limit_dim_size,
             var_list=expanded_var_list,
             empty_var_list=empty_var_list,  
            )
        
        # Finish up if no reprojection required
        dest_srs = get_spatial_ref_from_wkt(to_crs)
        if not to_crs or dest_srs.IsSame(get_spatial_ref_from_wkt(self.wkt)):
            logger.debug('No reprojection required for dataset {}'.format(nc_out_path))
            return
        
        try:
            logger.debug('Re-opening new dataset {}'.format(nc_out_path))
            new_dataset = netCDF4.Dataset(nc_out_path, 'r+')
            new_ncpu = NetCDFPointUtils(new_dataset, debug=self.debug)
            logger.debug('Reprojecting {} coordinates in new dataset'.format(len(new_ncpu.x_variable)))
            #TODO: Check coordinate variable data type if changing between degrees & metres
            new_ncpu._xycoords = transform_coords(new_ncpu.xycoords, self.wkt, to_crs)
            new_ncpu.x_variable[:] = new_ncpu._xycoords[:,0]
            new_ncpu.y_variable[:] = new_ncpu._xycoords[:,1]
            
            crs_variable_name, crs_variable_attributes = self.get_crs_attributes(to_crs)
            logger.debug('Setting {} variable attributes'.format(crs_variable_name))
            # Delete existing crs variable attributes
            for key in new_ncpu.crs_variable.__dict__.keys():
                if not key.startswith('_'):
                    delattr(new_ncpu.crs_variable, key) 
                    try:
                        delattr(new_ncpu.x_variable, key) 
                        delattr(new_ncpu.y_variable, key) 
                    except:
                        pass
                    
            # Set new crs variable attributes
            new_ncpu.crs_variable.setncatts(crs_variable_attributes)
            new_ncpu.x_variable.setncatts(crs_variable_attributes)
            new_ncpu.y_variable.setncatts(crs_variable_attributes)
            
            # Rename variables if switching between projected & unprojected
            if crs_variable_name != new_ncpu.crs_variable.name:
                logger.debug('Renaming {} variable to {}'.format(new_ncpu.crs_variable.name, crs_variable_name))
                new_dataset.renameVariable(new_ncpu.crs_variable.name, crs_variable_name) 
                
                if crs_variable_name == 'crs': # Geodetic
                    xy_varnames = ('longitude', 'latitude')
                    units = dest_srs.GetAngularUnitsName() + 's' # degrees
                    
                elif crs_variable_name in ['transverse_mercator', "albers_conical_equal_area"]: # Projected
                    xy_varnames = ('x', 'y')
                    units = units = dest_srs.GetLinearUnitsName() + 's' # metres
                    
                else:
                    raise BaseException('Unhandled crs variable name "{}"'.format(crs_variable_name))
                
                logger.debug('Renaming {} & {} variables to {} & {}'.format(new_ncpu.x_variable.name, 
                                                                             new_ncpu.y_variable.name,
                                                                             *xy_varnames
                                                                             ))
                new_dataset.renameVariable(new_ncpu.x_variable.name, xy_varnames[0])
                new_ncpu.x_variable.units = units
                new_ncpu.x_variable.long_name = xy_varnames[0]
                
                new_dataset.renameVariable(new_ncpu.y_variable.name, xy_varnames[1])
                new_ncpu.y_variable.units = units
                new_ncpu.y_variable.long_name = xy_varnames[1]
        
        finally:
            new_dataset.close()

    def set_global_attributes(self, compute_shape=False):
        '''\
        Function to set  global geometric metadata attributes in netCDF file
        N.B: This will fail if dataset is not writable
        '''
        try:
            metadata_srs = get_spatial_ref_from_wkt(METADATA_CRS)
            assert metadata_srs.IsGeographic(), 'Unable to set geodetic parameters for this dataset'
            
            
            #===================================================================
            # # Reopen as writable dataset
            # filepath = self.netcdf_dataset.filepath()
            # self.netcdf_dataset.close()
            # self.netcdf_dataset = netCDF4.Dataset(filepath, 'r+')
            #===================================================================
            
            logger.debug('Setting global geometric metadata attributes in netCDF point dataset with {} points'.format(self.netcdf_dataset.dimensions['point'].size))
            
            attribute_dict = dict(zip(['geospatial_lon_min', 'geospatial_lat_min', 'geospatial_lon_max', 'geospatial_lat_max'],
                              get_reprojected_bounds(self.bounds, self.wkt, METADATA_CRS)
                              )
                          )
            
            attribute_dict['geospatial_lon_units'] = 'degrees_east'
            attribute_dict['geospatial_lat_units'] = 'degrees_north'
            
            attribute_dict['geospatial_bounds_crs'] = metadata_srs.ExportToPrettyWkt()
            
            if compute_shape:
                try:
                    logger.debug('Computing concave hull')
                    attribute_dict['geospatial_bounds'] = shapely.wkt.dumps(
                        self.get_concave_hull(
                            to_wkt=METADATA_CRS,
                            buffer_distance=SHAPE_BUFFER_DISTANCE,
                            offset=SHAPE_OFFSET,
                            tolerance=SHAPE_SIMPLIFY_TOLERANCE,
                            max_polygons=SHAPE_MAX_POLYGONS,
                            max_vertices=SHAPE_MAX_VERTICES

                            ), 
                        rounding_precision=SHAPE_ORDINATE_DECIMAL_PLACES)
                except Exception as e:
                    logger.warning('Unable to compute concave hull shape: {}'.format(e))
                    try:
                        self.netcdf_dataset.geospatial_bounds = shapely.wkt.dumps(asPolygon([
                            [attribute_dict['geospatial_lon_min'], attribute_dict['geospatial_lat_min']], 
                            [attribute_dict['geospatial_lon_max'], attribute_dict['geospatial_lat_min']], 
                            [attribute_dict['geospatial_lon_max'], attribute_dict['geospatial_lat_max']],
                            [attribute_dict['geospatial_lon_min'], attribute_dict['geospatial_lat_max']], 
                            [attribute_dict['geospatial_lon_min'], attribute_dict['geospatial_lat_min']], 
                            ]))
                    except:
                        pass
                            
            logger.debug('attribute_dict = {}'.format(pformat(attribute_dict)))
        
            logger.debug('Writing global attributes to netCDF file'.format(self.netcdf_dataset.filepath()))
            for key, value in attribute_dict.items():
                setattr(self.netcdf_dataset, key, value)
                
            logger.debug('Finished setting global geometric metadata attributes in netCDF point dataset')
        except:
            logger.error('Unable to set geometric metadata attributes in netCDF point dataset')
            raise


    def set_variable_actual_range_attribute(self):
        '''\
        Function to set ACDD actual_range attribute in all non-index point-dimensioned variables
        N.B: Will fail if dataset is not writable
        '''
        self.netcdf_dataset.set_auto_mask(True)
        
        try:
            for variable_name, variable in self.netcdf_dataset.variables.items():
                # Skip all variables not of point dimensionality
                if 'point' not in variable.dimensions:
                    continue
                
                # Skip index variables
                if re.search('_index$', variable_name):
                    continue
                
                try:
                    variable.actual_range = np.array(
                        [np.nanmin(variable[:]), np.nanmax(variable[:])], dtype=variable.dtype)
                    logger.debug('{}.actual_range = {}'.format(variable_name, variable.actual_range))
                except:
                    logger.warning('Unable to compute actual_range value for point variable {}'.format(variable_name))
        except:
            logger.error('Unable to set variable actual_range metadata attributes in netCDF point dataset')
            raise



def main(debug=True):
    '''
    Main function for quick and dirty testing
    '''
    netcdf_path = sys.argv[1]
    
    netcdf_dataset = netCDF4.Dataset(netcdf_path, 'r')

    ncpu = NetCDFPointUtils(netcdf_dataset, debug=debug) # Enable debug output here
    
    # Create mask for last ten points
    mask = np.zeros(shape=(ncpu.point_count,), dtype='bool')
    mask[-10:] = True
    
    # Set list of fields to read
    field_list = None
    #field_list = ['latitude', 'longitude', 'obsno', 'reliab'] 
    
    point_data_generator = ncpu.all_point_data_generator(field_list, mask)
    
    # Retrieve point variable attributes first
    variable_attributes = next(point_data_generator)
    logger.info('variable_attributes: {}'.format(variable_attributes))

    # Use long names instead of variable names where they exist
    field_names = [variable_attributes[variable_name].get('long_name') or variable_name
                   for variable_name in variable_attributes.keys()]    
    logger.info('field_names: {}'.format(field_names))
    
    for point_data in point_data_generator:
        #logger.debug('point_data: {}'.format(pformat(point_data)))
        result_dict = dict(zip(field_names, point_data))
        logger.info('result_dict: {}'.format(result_dict))
        
       
if __name__ == '__main__':
    # Setup logging handlers if required
    console_handler = logging.StreamHandler(sys.stdout)
    #console_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)
    
    if not logger.handlers:
        # Set handler for root logger to standard output
        logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for logger {}'.format(logger.name))

    ncu_logger = logging.getLogger('geophys_utils._netcdf_utils')
    if not ncu_logger.handlers:
        ncu_logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for {}'.format(ncu_logger.name))

    main()        
