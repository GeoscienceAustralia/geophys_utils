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
Created on 14Sep.,2016

@author: Alex
'''
import numpy as np
import math
from scipy.ndimage import map_coordinates
from geophys_utils._crs_utils import get_utm_wkt, transform_coords
from geophys_utils._transect_utils import sample_transect
from geophys_utils._polygon_utils import netcdf2convex_hull
from geophys_utils._netcdf_utils import NetCDFUtils
from shapely.geometry import Polygon, MultiPolygon, asPolygon
from affine import Affine
import logging
import argparse
from distutils.util import strtobool
from shapely.geometry.base import BaseGeometry
import netCDF4
import sys
from skimage import measure

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module


class NetCDFGridUtils(NetCDFUtils):
    '''
    NetCDFGridUtils class to do various fiddly things with gridded NetCDF geophysics files.
    '''
    # Assume WGS84 lat/lon if no CRS is provided
    DEFAULT_CRS = "GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]]"
    HORIZONTAL_VARIABLE_NAMES = ['lon', 'Easting', 'x', 'longitude']
    DEFAULT_MAX_BYTES = 500000000  # Default to 500,000,000 bytes for NCI's OPeNDAP
    FLOAT_TOLERANCE = 0.000001

    def __init__(self, netcdf_dataset, debug=False):
        '''
        NetCDFGridUtils Constructor - wraps a NetCDF dataset
        '''
        def set_nominal_pixel_sizes():
            '''
            Function to set tuples with the nominal vertical and horizontal sizes of the centre pixel in metres and degrees
            '''
            centre_pixel_indices = [
                len(self.dimension_arrays[dim_index]) // 2 for dim_index in range(2)]

            # Get coordinates of centre pixel and next diagonal pixel
            centre_pixel_coords = [[self.dimension_arrays[dim_index][centre_pixel_indices[dim_index]] 
                                    for dim_index in range(2)],
                                   [self.dimension_arrays[dim_index][centre_pixel_indices[dim_index] + 1] 
                                    for dim_index in range(2)]
                                   ]

            if self.YX_order:
                for coord_index in range(2):
                    centre_pixel_coords[coord_index].reverse()

            #TODO: Make sure this is general for all CRSs
            self.x_variable = (self.netcdf_dataset.variables.get('longitude') 
                               or self.netcdf_dataset.variables.get('lon') 
                               or self.netcdf_dataset.variables.get('x')
                               )
            
            self.y_variable = (self.netcdf_dataset.variables.get('latitude') 
                               or self.netcdf_dataset.variables.get('lat') 
                               or self.netcdf_dataset.variables.get('y')
                               )
            
            self.y_inverted = (self.y_variable[-1] < self.y_variable[0]).item() # Should not have to deal with null values
        
            nominal_utm_wkt = get_utm_wkt(centre_pixel_coords[0], self.wkt)
            centre_pixel_utm_coords = transform_coords(
                centre_pixel_coords, from_wkt=self.wkt, to_wkt=nominal_utm_wkt)          
            
            self.nominal_pixel_metres = [round(abs(centre_pixel_utm_coords[1][
                        dim_index] - centre_pixel_utm_coords[0][dim_index]), 8) for dim_index in range(2)]
            
            
            centre_pixel_wgs84_coords = transform_coords(
                centre_pixel_coords, from_wkt=self.wkt, to_wkt='EPSG:4326')
            
            self.nominal_pixel_degrees = [round(abs(centre_pixel_wgs84_coords[1][
                        dim_index] - centre_pixel_wgs84_coords[0][dim_index]), 8) for dim_index in range(2)]

        def get_default_sample_metres():
            '''
            Function to return average nominal pixel size in metres rounded up to nearest 10^x or 5*10^x
            This is to provide a sensible default resolution for the sampling points along a transect by keeping it around the nominal pixel size
            '''
            log_10_avg_pixel_metres = math.log((self.nominal_pixel_metres[
                                               0] + self.nominal_pixel_metres[1]) / 2.0) / math.log(10.0)
            log_10_5 = math.log(5.0) / math.log(10.0)

            return round(math.pow(10.0, math.floor(log_10_avg_pixel_metres) +
                                  (log_10_5 if((log_10_avg_pixel_metres % 1.0) < log_10_5) else 1.0)))

    
        # Start of init function - Call inherited constructor first
        super().__init__(netcdf_dataset, debug=debug)
        
        logger.debug('Running NetCDFGridUtils constructor')
        
        self._GeoTransform = None


# assert len(self.netcdf_dataset.dimensions) == 2, 'NetCDF dataset must be
# 2D' # This is not valid

        try:
            data_variable_dimensions = [variable for variable in self.netcdf_dataset.variables.values() 
                                       if hasattr(variable, 'grid_mapping')][0].dimensions
            self._data_variable_list = [variable for variable in self.netcdf_dataset.variables.values() 
                                       if variable.dimensions == data_variable_dimensions]
        except:
            logger.debug('Unable to determine data variable(s) (must have same dimensions as variable with "grid_mapping" attribute)')
            raise
            
        #TODO: Make this work for multi-variate grids
        assert len(self.data_variable_list) > 0, 'Unable to determine any data variables (must have same dimensions as variable with "grid_mapping" attribute)'
        self.data_variable = self.data_variable_list[0]
        
        # Boolean flag indicating YX array ordering
        # TODO: Find a nicer way of dealing with this
        self.YX_order = self.data_variable.dimensions[
            1] in NetCDFGridUtils.HORIZONTAL_VARIABLE_NAMES

        # Two-element list of dimension varibles.
        self.dimension_arrays = [self.netcdf_dataset.variables[dimension_name][
            :] for dimension_name in self.data_variable.dimensions]

        self.pixel_size = [abs(self.GeoTransform[1]),
                           abs(self.GeoTransform[5])]
        
        self.pixel_count = list(self.data_variable.shape)
        
        if self.YX_order:
            self.pixel_size.reverse()
            self.pixel_count.reverse()

        self.min_extent = tuple([min(self.dimension_arrays[
                                dim_index]) - self.pixel_size[dim_index] / 2.0 for dim_index in range(2)])
        self.max_extent = tuple([max(self.dimension_arrays[
                                dim_index]) + self.pixel_size[dim_index] / 2.0 for dim_index in range(2)])

        set_nominal_pixel_sizes()

        self.default_sample_metres = get_default_sample_metres()

        # Create nested list of bounding box corner coordinates
        self.native_bbox = [[self.GeoTransform[0] + (x_pixel_offset * self.GeoTransform[1]) + (y_pixel_offset * self.GeoTransform[2]),
                             self.GeoTransform[3] + (x_pixel_offset * self.GeoTransform[4]) + (y_pixel_offset * self.GeoTransform[5])]
                            for x_pixel_offset, y_pixel_offset in [[0, self.pixel_count[1]], 
                                                                   [self.pixel_count[0], self.pixel_count[1]],
                                                                   [self.pixel_count[0], 0],
                                                                   [0, 0]
                                                                   ]
                            ]
        
        # Create bounds
        self.bounds = self.native_bbox[0] + self.native_bbox[2]

    def get_indices_from_coords(self, coordinates, wkt=None):
        '''
        Returns list of netCDF array indices corresponding to coordinates to support nearest neighbour queries
        @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
        @parameter wkt: Coordinate Reference System for coordinates. None == native NetCDF CRS
        '''
        wkt = wkt or self.wkt
        native_coordinates = transform_coords(coordinates, self.wkt, wkt)
        # Reshape 1D array into 2D single coordinate array if only one coordinate provided
        if native_coordinates.shape == (2,):
            native_coordinates = native_coordinates.reshape((1,2))
        
        # Convert coordinates to same dimension ordering as array
        if self.YX_order:
            native_coordinates = native_coordinates[:,1::-1]
                
        try:  # Multiple coordinates
            indices = [[np.where(abs(self.dimension_arrays[dim_index] - coordinate[dim_index]) <= (self.pixel_size[dim_index] / 2.0))[0][0] for dim_index in range(2)]
                       if not ([True for dim_index in range(2) if coordinate[dim_index] < self.min_extent[dim_index] or coordinate[dim_index] > self.max_extent[dim_index]])
                       else None
                       for coordinate in native_coordinates]
        except TypeError:  # Single coordinate pair
            indices = ([np.where(abs(self.dimension_arrays[dim_index] - native_coordinates[dim_index]) <= (self.pixel_size[dim_index] / 2.0))[0][0] for dim_index in range(2)]
                       if not [True for dim_index in range(2) if native_coordinates[dim_index] < self.min_extent[dim_index] or native_coordinates[dim_index] > self.max_extent[dim_index]]
                       else None)

        return indices

    def get_fractional_indices_from_coords(self, coordinates, wkt=None):
        '''
        Returns list of fractional array indices corresponding to coordinates to support interpolation
        @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
        @parameter wkt: Coordinate Reference System for coordinates. None == native NetCDF CRS
        '''
        wkt = wkt or self.wkt
        native_coordinates = transform_coords(coordinates, self.wkt, wkt)

        self.pixel_size

        # Convert coordinates to same order as array
        if self.YX_order:
            try:
                for coord_index in range(len(native_coordinates)):
                    if native_coordinates[coord_index] is not None:
                        native_coordinates[coord_index] = list(
                            native_coordinates[coord_index])
                        native_coordinates[coord_index].reverse()
            except:
                native_coordinates = list(native_coordinates)
                native_coordinates.reverse()
        # TODO: Make sure this still works with Southwards-positive datasets
        try:  # Multiple coordinates
            fractional_indices = [[(coordinate[dim_index] - min(self.dimension_arrays[dim_index])) / self.pixel_size[dim_index] for dim_index in range(2)]
                                  if not ([True for dim_index in range(2) if coordinate[dim_index] < self.min_extent[dim_index] or coordinate[dim_index] > self.max_extent[dim_index]])
                                  else None
                                  for coordinate in native_coordinates]
        except:  # Single coordinate pair
            fractional_indices = ([(native_coordinates[dim_index] - min(self.dimension_arrays[dim_index])) / self.pixel_size[dim_index] for dim_index in range(2)]
                                  if not [True for dim_index in range(2) if native_coordinates[dim_index] < self.min_extent[dim_index] or native_coordinates[dim_index] > self.max_extent[dim_index]]
                                  else None)

        return fractional_indices

    def get_value_at_coords(self, coordinates, wkt=None,
                            max_bytes=None, variable_name=None):
        '''
        Returns list of array values at specified coordinates
        @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
        @parameter wkt: WKT for coordinate Coordinate Reference System. None == native NetCDF CRS
        @parameter max_bytes: Maximum number of bytes to read in a single query. Defaults to NetCDFGridUtils.DEFAULT_MAX_BYTES
        @parameter variable_name: NetCDF variable_name if not default data variable
        '''
        # Use arbitrary maximum request size of NetCDFGridUtils.DEFAULT_MAX_BYTES
        # (500,000,000 bytes => 11180 points per query)
        #TODO: Find a better way of overcoming the netCDF problem where whole rows & columns are retrieved
        max_bytes = max_bytes or 100  # NetCDFGridUtils.DEFAULT_MAX_BYTES

        if variable_name:
            data_variable = self.netcdf_dataset.variables[variable_name]
        else:
            data_variable = self.data_variable

        no_data_value = data_variable._FillValue

        indices = np.array(self.get_indices_from_coords(coordinates, wkt))
        
#        return data_variable[indices[:,0], indices[:,1]].diagonal() # This could get too big

        # Allow for the fact that the NetCDF advanced indexing will pull back
        # n^2 cells rather than n
        max_points = max(
            int(math.sqrt(max_bytes / data_variable.dtype.itemsize)), 1)
        try:
            # Make this a vectorised operation for speed (one query for as many
            # points as possible)
            # Array of valid index pairs only
            index_array = np.array(
                [index_pair for index_pair in indices if index_pair is not None])
            assert len(index_array.shape) == 2 and index_array.shape[
                1] == 2, 'Not an iterable containing index pairs'
            # Boolean mask indicating which index pairs are valid
            mask_array = np.array([(index_pair is not None)
                                   for index_pair in indices])
            # Array of values read from variable
            value_array = np.ones(shape=(len(index_array)),
                                  dtype=data_variable.dtype) * no_data_value
            # Final result array including no-data for invalid index pairs
            result_array = np.ones(
                shape=(len(mask_array)), dtype=data_variable.dtype) * no_data_value
            start_index = 0
            end_index = min(max_points, len(index_array))
            while True:
                # N.B: ".diagonal()" is required because NetCDF doesn't do advanced indexing exactly like numpy
                # Hack is required to take values from leading diagonal. Requires n^2 elements retrieved instead of n. Not good, but better than whole array
                # TODO: Think of a better way of doing this
                value_array[start_index:end_index] = data_variable[
                    (index_array[start_index:end_index, 0], index_array[start_index:end_index, 1])].diagonal()
                if end_index == len(index_array):  # Finished
                    break
                start_index = end_index
                end_index = min(start_index + max_points, len(index_array))

            result_array[mask_array] = value_array
            return list(result_array)
        except:
            return data_variable[indices[0], indices[1]]

    def get_interpolated_value_at_coords(
            self, coordinates, wkt=None, max_bytes=None, variable_name=None):
        '''
        Returns list of interpolated array values at specified coordinates
        @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
        @parameter wkt: Coordinate Reference System for coordinates. None == native NetCDF CRS
        @parameter max_bytes: Maximum number of bytes to read in a single query. Defaults to NetCDFGridUtils.DEFAULT_MAX_BYTES
        @parameter variable_name: NetCDF variable_name if not default data variable
        '''
        # TODO: Check behaviour of scipy.ndimage.map_coordinates adjacent to no-data areas. Should not interpolate no-data value
        # TODO: Make this work for arrays > memory
        max_bytes = max_bytes or 100
        NetCDFGridUtils.DEFAULT_MAX_BYTES

        if variable_name:
            data_variable = self.netcdf_dataset.variables[variable_name]
        else:
            data_variable = self.data_variable

        no_data_value = data_variable._FillValue

        fractional_indices = self.get_fractional_indices_from_coords(
            coordinates, wkt)

        # Make this a vectorised operation for speed (one query for as many
        # points as possible)
        try:
            # Array of valid index pairs only
            index_array = np.array(
                [index_pair for index_pair in fractional_indices if index_pair is not None])
            assert len(index_array.shape) == 2 and index_array.shape[
                1] == 2, 'Not an iterable containing index pairs'
            # Boolean mask indicating which index pairs are valid
            mask_array = np.array([(index_pair is not None)
                                   for index_pair in fractional_indices])
            # Array of values read from variable
            value_array = np.ones(shape=(len(index_array)),
                                  dtype=data_variable.dtype) * no_data_value
            # Final result array including no-data for invalid index pairs
            result_array = np.ones(
                shape=(len(mask_array)), dtype=data_variable.dtype) * no_data_value

            value_array = map_coordinates(
                data_variable, index_array.transpose(), cval=no_data_value)

            result_array[mask_array] = value_array

            # Mask out any coordinates falling in no-data areas. Need to do this to stop no-data value from being interpolated
            # This is a bit ugly.
            result_array[np.array(self.get_value_at_coords(
                coordinates, wkt, max_bytes, variable_name)) == no_data_value] = no_data_value

            return list(result_array)
        except AssertionError:
            return map_coordinates(data_variable, np.array(
                [[fractional_indices[0]], [fractional_indices[1]]]), cval=no_data_value)


    def sample_transect(self, transect_vertices, wkt=None, sample_metres=None):
        '''
        Function to return a list of sample points sample_metres apart along lines between transect vertices
        @param transect_vertices: list or array of transect vertex coordinates
        @param wkt: coordinate reference system for transect_vertices
        @param sample_metres: distance between sample points in metres
        '''
        wkt = wkt or self.wkt
        sample_metres = sample_metres or self.default_sample_metres
        return sample_transect(transect_vertices, wkt, sample_metres)
        

    def get_convex_hull(self, to_wkt=None):
        '''\
        Function to return n x 2 array of coordinates for convex hull based on line start/end points
        Implements abstract base function in NetCDFUtils 
        @param to_wkt: CRS WKT for shape
        '''
        try:
            convex_hull = netcdf2convex_hull(self.netcdf_dataset, NetCDFGridUtils.DEFAULT_MAX_BYTES)
        except:
            logger.warning('Unable to compute convex hull. Using rectangular bounding box instead.')
            convex_hull = self.native_bbox
            
        return transform_coords(convex_hull, self.wkt, to_wkt)
    
    def get_concave_hull(self, 
                         to_wkt=None, 
                         buffer_distance=None, 
                         offset=None, 
                         tolerance=None, 
                         cap_style=1, 
                         join_style=1, 
                         max_polygons=10, 
                         max_vertices=1000
                         ):
        """\
        Returns the concave hull (as a shapely polygon) of grid edge points with data. 
        Implements abstract base function in NetCDFUtils 
        Note that all distance parameters are in pixel units, not in destination CRS units
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
        PAD_WIDTH = 1
        
        # Parameters for shape simplification in pixel sizes
        buffer_distance = buffer_distance or max(self.data_variable.shape) // 20 # Tune this to suit overall size of grid
        offset = offset or 0.5 # Take final shape half a pixel out from centre coordinates to cover pixel edges
        tolerance = tolerance or 0.5        
        
        
        def discard_internal_polygons(geometry):
            '''\
            Helper function to discard internal polygons
            '''
            if type(geometry) == MultiPolygon:
                polygon_list = []
                for polygon in geometry:
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
                    external_geometry = polygon_list[0] # Single polygon
                else:
                    external_geometry = MultiPolygon(polygon_list)
                    
                logger.debug('{} internal polygons discarded.'.format(len(geometry) - len(polygon_list)))
                    
            elif type(geometry) == Polygon:
                external_geometry = Polygon(geometry.exterior)
            else:
                raise ValueError('Unexpected type of geometry: {}'.format(type(geometry)))
            
            return external_geometry
            
        def contour_to_polygon(polygon_vertices):
            '''\
            Helper function to turn an individual contour into shapely Polygon in the correct xy axis order
            @param polygon_vertices: n x 2 array of pixel coordinates forming a polygon
            '''
            # Note that we need coordinates in xy order, but the polygon vertices are in array order (probably yx)
            if self.YX_order: # yx array order - y-axis flip required
                reorder_slices = (slice(None, None, None), slice(None, None, -1))
            else: # xy array order - no y-axis flip required
                reorder_slices = (slice(None, None, None), slice(None, None, None))
                
            # Create polygon in pixel ordinates, and then buffer outwards by 0.5 pixel widths, and simplify by 0.25 pixel widths
            # Note that this polygon is in array order, i.e. probably yx and not xy, so we need to apply the reorder_slices
            return asPolygon(polygon_vertices[reorder_slices])
        
        def transform_geometry_pixel_to_wkt(geometry, to_wkt):
            '''\
            Helper function to transform geometry from pixel coordinates to specified WKT via (self.wkt)
            # N.B: Ignores polygon interiors
            @param geometry: Shapely MultiPolygon or Polygon to transform
            @param to_wkt: WKT of destination CRS
            '''
            # Set affine transform from GeoTransform values for CRS given by self.wkt
            affine_transform = Affine.from_gdal(*self.GeoTransform)
        
            if type(geometry) == MultiPolygon:
                return MultiPolygon([
                    asPolygon(transform_coords(np.array([(affine_transform * pixel_coordinate_pair) 
                                                    for pixel_coordinate_pair in polygon.exterior.coords
                                                    ]),
                         self.wkt, to_wkt)
                         )
                    for polygon in geometry
                    ])

            elif type(geometry) == Polygon:
                return asPolygon(transform_coords(np.array([(affine_transform * pixel_coordinate_pair) 
                                                    for pixel_coordinate_pair in geometry.exterior.coords
                                                    ]),
                         self.wkt, to_wkt)
                         )
            else:
                raise ValueError('Unexpected type of geometry: {}'.format(type(geometry)))


        def get_offset_geometry(geometry, buffer_distance, offset, tolerance, cap_style, join_style, max_polygons, max_vertices):
            '''\
            Helper function to return offset geometry. Will keep trying larger buffer_distance values until there is a manageable number of polygons
            '''
            logger.debug('Computing offset geometry with buffer_distance = {}'.format(buffer_distance))
            offset_geometry = (geometry.buffer(
                buffer_distance, 
                cap_style=cap_style, 
                join_style=join_style).simplify(tolerance).buffer(
                    offset-buffer_distance, 
                    cap_style=cap_style, 
                    join_style=join_style).simplify(tolerance)
                )      
            discard_internal_polygons(offset_geometry)
            
            # Keep doubling the buffer distance if there are too many polygons
            if (
                (max_polygons and type(offset_geometry) == MultiPolygon and len(offset_geometry) > max_polygons)
                or
                (max_vertices and type(offset_geometry) == MultiPolygon and 
                    sum([len(polygon.exterior.coords) #+ sum([len(interior_ring.coords) for interior_ring in polygon.interiors]) 
                         for polygon in offset_geometry]) > max_vertices)
                or
                (max_vertices and type(offset_geometry) == Polygon and 
                    (len(offset_geometry.exterior.coords) #+ sum([len(interior_ring.coords) for interior_ring in offset_geometry.interiors])
                     ) > max_vertices)
                ):
                return get_offset_geometry(geometry, buffer_distance*2, offset, tolerance, cap_style, join_style, max_polygons, max_vertices)
                
            return offset_geometry
        # Create data/no-data mask
        mask = (self.data_variable[:] != self.data_variable._FillValue)

        # pad with nodata so that boundary edges are detected
        padded_mask = np.pad(mask, pad_width=PAD_WIDTH, mode='constant', constant_values=False)

        # find contours where the high pieces (data) are fully connected
        # that there are no unnecessary holes in the polygons
        # shift the coordinates back by 1 to get original unpadded pixel coordinates
        logger.debug('Generating contours for grid of size {}'.format(', '.join(str(size) for size in self.data_variable.shape)))
        contours = [contour - np.array([[PAD_WIDTH, PAD_WIDTH]])
                    for contour in measure.find_contours(padded_mask, 0.5, fully_connected='high')]
        #logger.debug('contours = {}'.format(contours))
        
        concave_hull = MultiPolygon([contour_to_polygon(contour) for contour in contours]).simplify(tolerance)
        concave_hull = get_offset_geometry(concave_hull, buffer_distance, offset, tolerance, cap_style, join_style, max_polygons, max_vertices)
        
        if self.debug:
            if type(concave_hull) == MultiPolygon:
                polygon_count = len(concave_hull)
                vertex_count = sum([len(polygon.exterior.coords) for polygon in concave_hull])

            elif type(concave_hull) == Polygon:
                polygon_count = 1
                vertex_count = len(concave_hull.exterior.coords)
                
            else:
                raise ValueError('Unexpected type of geometry: {}'.format(type(concave_hull)))
            logger.debug('Final shape has {} vertices in {} polygons.'.format(vertex_count, polygon_count))
            

        # Transform from pixel ordinates to to_wkt (via self.wkt using affine_transform)
        return transform_geometry_pixel_to_wkt(concave_hull, to_wkt)  
    
    def get_dimension_ranges(self, bounds, bounds_wkt=None):
        '''
        Function to dict of (start, end+1) tuples keyed by dimension name from a bounds geometry or ordinates
        @parameter bounds: Either an iterable containing [<xmin>, <ymin>, <xmax>, <ymax>] or a shapely (multi)polygon
        @parameter bounds_wkt: WKT for bounds CRS. Defaults to dataset native CRS
        @return dim_range_dict: dict of (start, end+1) tuples keyed by dimension name
        '''
        if isinstance(bounds, BaseGeometry): # Process shapely (multi)polygon bounds        
            bounds_ordinates = bounds.bounds # Obtain [<xmin>, <ymin>, <xmax>, <ymax>] from Shapely geometry
        else:
            bounds_ordinates = bounds # Use provided [<xmin>, <ymin>, <xmax>, <ymax>] parameter
            
        if bounds_wkt is not None: # Reproject bounds to native CRS if required
            bounds_ordinates = self.get_reprojected_bounds(bounds_ordinates, bounds_wkt, self.wkt)

        #=======================================================================
        # bounds_half_size = abs(np.array([bounds_ordinates[2] - bounds_ordinates[0], bounds_ordinates[3] - bounds_ordinates[1]])) / 2.0
        # bounds_centroid = np.array([bounds_ordinates[0], bounds_ordinates[1]]) + bounds_half_size
        #=======================================================================
        
        dimension_names = self.data_variable.dimensions
            
        dim_range_dict = {}
        try:
            logger.debug('self.dimension_arrays = {}'.format(self.dimension_arrays))
            for dim_index in range(2):
                #TODO: Maybe make this work for pixel edges, not centres
                #TODO: Fix indexing in self.dimension_arrays where 1-dim_index only works for lon-lat array
                subset_indices = np.where(np.logical_and(self.dimension_arrays[1-dim_index] >= bounds_ordinates[dim_index],
                                                         self.dimension_arrays[1-dim_index] <= bounds_ordinates[dim_index+2]))[0]
                                                         
                logger.debug('subset_indices = {}'.format(subset_indices))
                
                logger.debug('self.dimension_arrays[{}].shape = {}'.format(1-dim_index, self.dimension_arrays[1-dim_index].shape))
                                                         
                dim_range_dict[dimension_names[1-dim_index]] = (subset_indices[0], min(subset_indices[-1]+1,
                                                                                     self.dimension_arrays[1-dim_index].shape[0]
                                                                                     ) # add 1 to upper index
                                                              )
                
                logger.debug('dim_range_dict["{}"] = {}'.format(dimension_names[1-dim_index], dim_range_dict[dimension_names[1-dim_index]]))
            
            return dim_range_dict
        except Exception as e:
            logger.debug('Unable to determine range indices: {}'.format(e))
            return None

    
    @property
    def GeoTransform(self):
        '''
        Property getter function to return geotransform as required
        '''
        if not self._GeoTransform:
            try:
                # Assume string or array representation of GeoTransform exists
                self._GeoTransform = self.crs_variable.GeoTransform
            except:
                #TODO: create GeoTransform from x & y variables
                raise BaseException('Unable to determine GeoTransform')

            if type(self._GeoTransform) == str:
                # Convert string representation of GeoTransform to array
                self._GeoTransform = [float(number.strip())
                                      for number in self.crs_variable.GeoTransform.strip().split(' ')
                                      ]
                   
        return self._GeoTransform

    def copy(self, 
             nc_out_path, 
             datatype_map_dict={},
             variable_options_dict={},
             dim_range_dict={},
             dim_mask_dict={},
             nc_format=None,
             limit_dim_size=False,
             empty_var_list=[],
             invert_y=None):
        '''
        Function to copy a netCDF dataset to another one with potential changes to size, format, 
            variable creation options and datatypes.
            
            @param nc_out_path: path to netCDF output file 
            @param datatype_map_dict: dict containing any maps from source datatype to new datatype.
                e.g. datatype_map_dict={'uint64': 'uint32'}  would convert all uint64 variables to uint32.
            @param variable_options_dict: dict containing any overrides for per-variable variable creation 
                options. e.g. variable_options_dict={'sst': {'complevel': 2, 'zlib': True}} would apply
                compression to variable 'sst'
            @param dim_range_dict: dict of (start, end+1) tuples keyed by dimension name
            @param dim_mask_dict: dict of boolean arrays keyed by dimension name
            @param nc_format: output netCDF format - 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET', 
                'NETCDF3_64BIT_DATA', 'NETCDF4_CLASSIC', or 'NETCDF4'. Defaults to same as input format.  
            @param limit_dim_size: Boolean flag indicating whether unlimited dimensions should be fixed
            @param empty_var_list: List of strings denoting variable names for variables which should be created but not copied
            @param invert_y: Boolean parameter indicating whether copied Y axis should be Southwards positive (None means same as source)
        '''  
        # Call inherited NetCDFUtils method
        super().copy( 
             nc_out_path, 
             datatype_map_dict=datatype_map_dict,
             variable_options_dict=variable_options_dict,
             dim_range_dict=dim_range_dict,
             dim_mask_dict=dim_mask_dict,
             nc_format=nc_format,
             limit_dim_size=limit_dim_size,
             empty_var_list=empty_var_list
             )
        try:
            logger.debug('Re-opening new dataset {}'.format(nc_out_path))
            new_dataset = netCDF4.Dataset(nc_out_path, 'r+')
            new_ncgu = NetCDFGridUtils(new_dataset, debug=self.debug)
            
            if dim_range_dict: # Subsets were requested
                logger.debug('dim_range_dict = {}'.format(dim_range_dict))
                logger.debug('Old first coordinate = ({}, {})'.format(self.x_variable[0], self.y_variable[0]))
                logger.debug('New first coordinate = ({}, {})'.format(new_ncgu.x_variable[0], new_ncgu.y_variable[0]))
                geo_transform = new_ncgu.GeoTransform
                logger.debug('Old GeoTransform = {}'.format(geo_transform))
                for gt_ordinate_index, ordinate_variable, gt_pixel_width_index in [(0, new_ncgu.x_variable, 1), (3, new_ncgu.y_variable, 5)]:
                    if ordinate_variable.dimensions[0] in dim_range_dict.keys():
                        geo_transform[gt_ordinate_index] = ordinate_variable[0].item() - (geo_transform[gt_pixel_width_index] / 2.0) # Subtract half pixel width
                logger.debug('New GeoTransform = {}'.format(geo_transform))
                gt_string = ' '.join([str(value) for value in geo_transform])
                logger.info('Updating GeoTransform to "{}"'.format(gt_string))
                new_ncgu.crs_variable.GeoTransform = gt_string
            
            if invert_y is None: # No change required
                return 
            
            if invert_y == new_ncgu.y_inverted:
                logger.debug('{} does not require any alteration to make y-axis inversion {}'.format(nc_out_path, invert_y))
                return
            
            assert new_ncgu.y_variable, 'No Y-axis indexing variable defined in {}'.format(nc_out_path)
            
            assert len(new_ncgu.y_variable.shape) == 1, 'Y-axis indexing variable should be 1D'
            
            y_axis_dimension_name = new_ncgu.y_variable.dimensions[0]
            
            for variable_name, variable in new_dataset.variables.items():
                if y_axis_dimension_name in variable.dimensions:
                    logger.debug('Inverting Y dimension of variable {}'.format(variable_name))
                    variable_slices = [
                        slice(None, None, -1) if dimension_name == y_axis_dimension_name else slice(None, None, None)
                        for dimension_name in variable.dimensions
                        ]
                    variable[:] = variable[variable_slices]
                    
        finally:
            new_dataset.close()
                       

    
def main():
    '''
    Main function for quick and dirty testing
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--format", help="NetCDF file format (one of 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA')",
                        type=str, default='NETCDF4')
    parser.add_argument("--chunkspec", help="comma-separated list of <dimension_name>/<chunk_size> specifications",
                        type=str)
    parser.add_argument("--complevel", help="Compression level for chunked variables as an integer 0-9. Default is 4",
                        type=int, default=4)
    parser.add_argument('-i', '--invert_y', help='Store copy with y-axis indexing Southward positive (grids only)', type=str)
    parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    
    args = parser.parse_args()
    
    if args.invert_y is not None:
        invert_y = bool(strtobool(args.invert_y))
    else:
        invert_y = None # Default to same as source
    
    if args.chunkspec:
        chunk_spec = {dim_name: int(chunk_size) 
                    for dim_name, chunk_size in [chunk_spec_string.strip().split('/') for chunk_spec_string in args.chunkspec.split(',')]}
    else:
        chunk_spec = None
            
    ncgu = NetCDFGridUtils(args.input_path,
                      debug=args.debug
                      )   
    
    ncgu.copy(args.output_path, 
             #datatype_map_dict={},
             # Compress all chunked variables
             variable_options_dict={variable_name: {'chunksizes': [chunk_spec.get(dimension) 
                                                                   for dimension in variable.dimensions
                                                                   ],
                                                    'zlib': bool(args.complevel),
                                                    'complevel': args.complevel
                                                    }
                               for variable_name, variable in ncgu.netcdf_dataset.variables.items()
                               if (set(variable.dimensions) & set(chunk_spec.keys()))
                               } if chunk_spec else {},
             #dim_range_dict={'lat': (5,205),'lon': (5,305)},
             #dim_mask_dict={},
             nc_format=args.format,
             #limit_dim_size=False
             invert_y=invert_y
             )
        
if __name__ == '__main__':
    console_handler = logging.StreamHandler(sys.stdout)
    # console_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(name)s: %(message)s')
    console_handler.setFormatter(console_formatter)
 
    if not logger.handlers:
        logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for {}'.format(logger.name))

    main()
