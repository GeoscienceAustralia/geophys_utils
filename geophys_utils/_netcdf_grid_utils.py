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
from geophys_utils._polygon_utils import netcdf2convex_hull, get_grid_edge_points
from geophys_utils._netcdf_utils import NetCDFUtils
from geophys_utils._concave_hull import concaveHull
from shapely.geometry import shape, Polygon, MultiPolygon, asMultiPoint
import logging
import argparse
from distutils.util import strtobool
from shapely.geometry.base import BaseGeometry
import gc

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
            self.y_variable = (self.netcdf_dataset.variables.get('lat') 
                               or self.netcdf_dataset.variables.get('y')
                               )
            
            self.y_inverted = (self.y_variable[-1] < self.y_variable[0])
        
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
        assert len(self.data_variable_list) == 1, 'Unable to determine single data variable (must have "grid_mapping" attribute)'
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
            #logger.info('Unable to compute convex hull. Using rectangular bounding box instead.')
            convex_hull = self.native_bbox
            
        return transform_coords(convex_hull, self.wkt, to_wkt)
    
    def get_concave_hull(self, 
                         to_wkt=None, 
                         buffer_distance=None, 
                         offset=0, 
                         tolerance=0.0005, 
                         cap_style=1, 
                         join_style=1, 
                         max_polygons=10, 
                         max_vertices=1000
                         ):
        """\
        Returns the concave hull (as a shapely polygon) of grid edge points with data. 
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
        def get_offset_geometry(geometry, buffer_distance, offset, tolerance, cap_style, join_style, max_polygons, max_vertices):
            '''\
            Helper function to return offset geometry. Will keep trying larger buffer_distance values until there is a manageable number of polygons
            '''
            logger.debug('Computing offset geometry with buffer_distance = {}'.format(buffer_distance))
            offset_geometry = geometry.buffer(buffer_distance, cap_style=cap_style, join_style=join_style).simplify(tolerance)
            offset_geometry = offset_geometry.buffer(offset-buffer_distance, cap_style=cap_style, join_style=join_style).simplify(tolerance)

            if type(offset_geometry) == MultiPolygon:
                offset_geometry = MultiPolygon([Polygon(p.exterior) for p in offset_geometry])
            elif type(offset_geometry) == Polygon:
                offset_geometry = Polygon(offset_geometry.exterior)
            else:
                raise ValueError('Unexpected type of geometry: {}'.format(type(offset_geometry)))
            
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

        buffer_distance = buffer_distance or 2.0 * max(*self.pixel_size) # Set initial buffer_distance to 2 x pixel size in native units
        logger.debug('Initial buffer_distance = {}'.format(buffer_distance))

        edge_multipoint = asMultiPoint(np.array(get_grid_edge_points(self.data_variable, self.dimension_arrays, self.data_variable._FillValue)))        
        
        offset_geometry = get_offset_geometry(edge_multipoint, buffer_distance, 0, tolerance, cap_style, join_style, max_polygons, max_vertices)
    
        # Transform vertices from native CRS to required CRS
        if type(offset_geometry) == MultiPolygon:
            offset_geometry = MultiPolygon([Polygon(transform_coords(polygon.exterior.coords, self.wkt, to_wkt)) for polygon in offset_geometry])
        elif type(offset_geometry) == Polygon:
            offset_geometry = Polygon(transform_coords(offset_geometry.exterior.coords, self.wkt, to_wkt))
        else:
            raise ValueError('Unexpected type of geometry: {}'.format(type(offset_geometry)))
        
        return offset_geometry

    
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


    
def main():
    '''
    Main function for quick and dirty testing
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-c', '--copy', 
                        dest='do_copy', 
                        action='store_const', 
                        const=True, default=False,
                        help='Copy netCDF files')
    parser.add_argument("-f", "--format", help="NetCDF file format (one of 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA')",
                        type=str, default='NETCDF4')
    parser.add_argument("--chunkspec", help="comma-separated list of <dimension_name>/<chunk_size> specifications",
                        type=str)
    parser.add_argument("--complevel", help="Compression level for chunked variables as an integer 0-9. Default is 4",
                        type=int, default=4)
    parser.add_argument('-i', '--invert_y', help='Store copy with y-axis indexing Southward positive', type=str)
    parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    
    args = parser.parse_args()
    
    if args.invert_y is not None:
        invert_y = bool(strtobool(args.invert_y))
    else:
        invert_y = None # Default to same as source
    
    if args.do_copy:
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
             #dim_range_dict={},
             nc_format=args.format,
             #limit_dim_size=False
             invert_y=invert_y
             )
        

if __name__ == '__main__':
    main()
