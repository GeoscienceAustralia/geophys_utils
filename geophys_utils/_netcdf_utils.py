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
NetCDFUtils class implementing useful functionality against netCDF files

Created on 2Mar.,2017

@author: u76345
'''
import abc
import netCDF4
import math
import itertools
import argparse
import re
import sys
import numpy as np
import logging
import osgeo
from pprint import pformat
from geophys_utils._crs_utils import transform_coords, get_spatial_ref_from_wkt

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

class NetCDFUtils(object):
    '''
    NetCDFUtils class implementing useful functionality against netCDF files
    '''
    # Point, line  and grid subclasses will need this, even though this class is non-spatial
    X_DIM_VARIABLE_NAMES = ['longitude', 'lon', 'x', 'Easting']
    Y_DIM_VARIABLE_NAMES = ['latitude', 'lat', 'y', 'Northing']
    CRS_VARIABLE_NAMES = ['crs', 'transverse_mercator']
    
    DEFAULT_COPY_OPTIONS = {'complevel': 4, 
                            'zlib': True, 
                            'fletcher32': True,
                            'shuffle': True,
                            'endian': 'little',
                            #'chunksizes': [1024, 1024]
                            }

    def __init__(self, netcdf_dataset, debug=False):
        '''
        Constructor for NetCDFUtils
        '''
        self._debug = None # Initialise private variable
        self.debug = debug # Set debug property
        
        logger.debug('Running NetCDFUtils constructor')
        
        if type(netcdf_dataset) == str: # String provided as path to netCDF file
            self.nc_path = netcdf_dataset
            self._netcdf_dataset = None
        elif type(netcdf_dataset) == netCDF4.Dataset: # NetCDF4.Dataset object provided
            self._netcdf_dataset = netcdf_dataset 
            self.nc_path = netcdf_dataset.filepath()
        else:
            raise TypeError('Invalid netcdf_dataset type')
        
        self.opendap = (re.match('^http.*', self.nc_path) is not None)
        if self.opendap:
            self.max_bytes = 500000000 # 500MB limit for NCI OPeNDAP
        else:
            self.max_bytes = 4000000000 # 4GB limit for direct netCDF file access
        
        # Initialise private property variables to None until set by property getter methods
        self._data_variable_list = None       
        self._crs_variable = None # Needs to be set in subclass constructor
        self._wkt = None
        self._wgs84_bbox = None
        
        self._x_variable = None
        for variable_name in NetCDFUtils.X_DIM_VARIABLE_NAMES:
            self._x_variable = self.netcdf_dataset.variables.get(variable_name) 
            if self._x_variable is not None:
                break
        
        self._y_variable = None
        for variable_name in NetCDFUtils.Y_DIM_VARIABLE_NAMES:
            self._y_variable = self.netcdf_dataset.variables.get(variable_name) 
            if self._y_variable is not None:
                break
        
        # Set auto masking to True for consistent behaviour
        self.netcdf_dataset.set_auto_mask(True) #TODO: set this at a function level
        

    def copy(self, 
             nc_out_path, 
             datatype_map_dict={},
             variable_options_dict={},
             dim_range_dict={},
             dim_mask_dict={},
             nc_format=None,
             limit_dim_size=False,
             empty_var_list=[]):
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
        '''  
        logger.debug('variable_options_dict: {}'.format(variable_options_dict))   
                  
        self.netcdf_dataset.set_auto_mask(False)
        
        # Override default variable options with supplied ones for all data variables
        for variable_name in self._netcdf_dataset.variables.keys():
            variable_dict = dict(NetCDFUtils.DEFAULT_COPY_OPTIONS)
            variable_dict.update(variable_options_dict.get(variable_name) or {})
            variable_options_dict[variable_name] = variable_dict
                                
        nc_format = nc_format or self.netcdf_dataset.file_format 
        logger.debug('Output format is %s' % nc_format)
        
        nc_output_dataset = netCDF4.Dataset(nc_out_path, mode="w", clobber=True, format=nc_format)
        nc_output_dataset.set_auto_mask(False)
        
        try:
            dims_used = set()
            dim_size = {}
            for variable_name, variable in self.netcdf_dataset.variables.items():
                dims_used |= set(variable.dimensions)
                
                # Update the sizes for all dimensions which have masks or ranges
                for dimension_index in range(len(variable.dimensions)): 
                    dimension_name = variable.dimensions[dimension_index]
                    source_dimension = self.netcdf_dataset.dimensions[dimension_name]
                    
                    if dim_size.get(dimension_name): # We have already configured this dimension
                        continue
                        
                    dim_mask = dim_mask_dict.get(dimension_name)
                    if dim_mask is None:
                        dim_mask = np.ones(shape=(variable.shape[dimension_index],), dtype=np.bool)
                    else:
                        assert dim_mask.shape == (source_dimension.size,), 'Dimension mask must be a 1D boolean mask of size {}'.format(source_dimension.size)
                        
                    dim_range = dim_range_dict.get(dimension_name)
                    if dim_range:
                        dim_mask[:dim_range[0]] = False
                        dim_mask[dim_range[1]:] = False
                        
                    dim_size[dimension_name] = np.count_nonzero(dim_mask) # Update sizes to take masks into account
                        
                    #dim_mask_dict[dimension_name] = dim_mask # Update mask to include range
                        
                    
            #logger.debug(dim_size)
            
            #Copy dimensions
            for dimension_name, dimension in self.netcdf_dataset.dimensions.items():
                if dimension_name in dims_used: # Discard unused dimensions
                    logger.debug('Copying dimension %s of length %d' % (dimension_name, dim_size[dimension_name]))
                    nc_output_dataset.createDimension(dimension_name, 
                                          dim_size[dimension_name] 
                                          if not dimension.isunlimited() or limit_dim_size 
                                          else None)
                else:
                    logger.debug('Skipping unused dimension %s' % dimension_name)
    
            # Copy variables
            for variable_name, input_variable in self.netcdf_dataset.variables.items():
                dtype = datatype_map_dict.get(str(input_variable.datatype)) or input_variable.datatype
                
                # Special case for "crs" or "transverse_mercator" - want byte datatype
                if input_variable == self.crs_variable: 
                    dtype = 'i1'
                    
                # Start off by copying options from input variable (if specified)
                var_options = input_variable.filters() or {}
                
                # Chunking is defined outside the filters() result
                chunking = input_variable.chunking()
                if chunking and chunking != 'contiguous':
                    # Input variable is chunked - use same chunking by default unless overridden
                    input_variable_chunking = [min(chunking[dimension_index], dim_size[input_variable.dimensions[dimension_index]])
                                                 for dimension_index in range(len(chunking))]
                elif (len(input_variable.dimensions) == 2 and 
                    variable_options_dict.get(variable_name) and
                    variable_options_dict.get(variable_name).get('chunksizes')
                    ): #TODO: Improve this
                    # If input variable is unchunked 2D and output chunking is specified - assume row chunking for input
                    input_variable_chunking = [1, dim_size[input_variable.dimensions[1]]]
                else:
                    # Input variable is not chunked
                    input_variable_chunking = None
                
                # Default to same chunking on input and output
                if input_variable_chunking: 
                    var_options['chunksizes'] = input_variable_chunking
                    
                if hasattr(input_variable, '_FillValue'):
                    var_options['fill_value'] = input_variable._FillValue
                    
                # Apply any supplied options over top of defaults
                var_options.update(variable_options_dict.get(variable_name) or {})
                
                # Ensure chunk sizes aren't bigger than variable sizes
                if var_options.get('chunksizes'):
                    for dimension_index in range(len(input_variable.dimensions)):
                        var_options['chunksizes'][dimension_index] = min(var_options['chunksizes'][dimension_index] or dim_size[input_variable.dimensions[dimension_index]],
                                                                         dim_size[input_variable.dimensions[dimension_index]])
                             
                options_string = ' with options: %s' % ', '.join(['%s=%s' % item for item in var_options.items()]) if var_options else ''   
                logger.debug("Copying variable %s from datatype %s to datatype %s%s" % (variable_name, 
                                                                                       input_variable.datatype, 
                                                                                       dtype, 
                                                                                       options_string
                                                                                       )
                            )
                # Create output variable using var_options to specify output options
                output_variable = nc_output_dataset.createVariable(variable_name, 
                                              dtype, 
                                              input_variable.dimensions,
                                              **var_options
                                              )
                
                # Copy variable attributes
                logger.debug('\tCopying %s attributes: %s' % (variable_name, ', '.join(input_variable.ncattrs())))
                output_variable.setncatts({k: input_variable.getncattr(k) for k in input_variable.ncattrs() if not k.startswith('_')})
                
                if variable_name not in empty_var_list:
                    # Copy data
                    if input_variable.shape: # array
                        input_variable_slices = tuple([
                            slice(*dim_range_dict[input_variable.dimensions[dimension_index]])  
                            if dim_range_dict.get(input_variable.dimensions[dimension_index])
                            else slice(0, input_variable.shape[dimension_index])
                        for dimension_index in range(len(input_variable.dimensions))
                        ])
                        
                        logger.debug('input_variable_slices={}'.format(input_variable_slices))
                        
                        logger.debug('\tCopying {} array data of shape {}'.format(
                            variable_name,
                            tuple([input_variable_slices[dimension_index].stop - input_variable_slices[dimension_index].start
                                   for dimension_index in range(len(input_variable.dimensions))]
                                  )
                            ))
                        
                        # Build list of full dimension masks for this variable
                        # We are using None instead of all-true masks because of the issue described in 
                        # https://github.com/numpy/numpy/issues/13255
                        # This results in extra dimensions, but the piece assignment operation doesn't care
                        # See the behaviour described in https://stackoverflow.com/questions/1408311/numpy-array-slice-using-none
                        variable_masks = tuple([dim_mask_dict.get(dimension_name) # Use specified mask if it exists, None otherwise
                                          for dimension_name in input_variable.dimensions
                                          ])
                        
                        if not input_variable_chunking: 
                            # No chunking - Try to copy in one hit. This may bork due to OPeNDAP or memory limitations
                            #TODO: Make this safe for massive arrays, possibly using the array_pieces code
                            if any(variable_masks):
                                output_variable[...] = input_variable[input_variable_slices][variable_masks]
                            else:
                                output_variable[...] = input_variable[input_variable_slices]
                        
                        else: # Chunked - perform copy in pieces
                            #TODO: Improve performance for small chunks, and maybe look at chunk alignment for slices
                             
                            # Use largest chunk sizes between input and output
                            piece_sizes = [
                                max(var_options['chunksizes'][dimension_index],
                                    input_variable_chunking[dimension_index]
                                    )
                                for dimension_index in range(len(input_variable.dimensions))
                                ]     
                                                   
                            
                            piece_index_ranges = [
                                (input_variable_slices[dimension_index].start // piece_sizes[dimension_index],
                                 int(math.ceil(float(input_variable_slices[dimension_index].stop) / piece_sizes[dimension_index]))
                                 )
                                for dimension_index in range(len(input_variable.dimensions))
                                ]
                                                  
                            piece_counts = [
                                int(math.ceil(float(dim_size[input_variable.dimensions[dimension_index]]) / piece_sizes[dimension_index]))
                                for dimension_index in range(len(input_variable.dimensions)) 
                                ]
                         
                            logger.debug('\tCopying {} pieces of dimensions {}'.format(' x '.join([str(piece_count) for piece_count in piece_counts]),
                                                                                       ' x '.join([str(piece_size) for piece_size in piece_sizes])
                                                                          )
                                        )
                             
                            # Iterate over every piece
                            for piece_indices in itertools.product(
                                *[range(piece_index_ranges[dimension_index][0], 
                                        piece_index_ranges[dimension_index][1]
                                        )
                                  for dimension_index in range(len(input_variable.dimensions))
                                  ]
                                 ):
                                 
                                offset_piece_indices = tuple([(indices[0] - indices[1] + 1) 
                                                              for indices in zip(piece_indices, [piece_index_range[0] for piece_index_range in piece_index_ranges])])
                                logger.debug('\t\tCopying piece {}'.format(offset_piece_indices))
                                
                                 
                                piece_read_slices = tuple([
                                    slice(
                                        max(input_variable_slices[dimension_index].start,
                                            piece_indices[dimension_index] * piece_sizes[dimension_index]
                                            ),                                                 
                                        min(input_variable_slices[dimension_index].stop,
                                            (piece_indices[dimension_index] + 1) * piece_sizes[dimension_index]
                                            )
                                        )
                                    for dimension_index in range(len(input_variable.dimensions))
                                    ])
                                 
                                logger.debug('piece_read_slices = {}'.format(piece_read_slices))
                                
                                piece_write_slices = tuple([
                                    (
                                        slice(
                                            np.count_nonzero(variable_masks[dimension_index][input_variable_slices[dimension_index]][:piece_read_slices[dimension_index].start]),
                                            (np.count_nonzero(variable_masks[dimension_index][input_variable_slices[dimension_index]][:piece_read_slices[dimension_index].start]) + 
                                                np.count_nonzero(variable_masks[dimension_index][input_variable_slices[dimension_index]][piece_read_slices[dimension_index]]))
                                            ) 
                                        if variable_masks[dimension_index] is not None else 
                                        slice(piece_read_slices[dimension_index].start - input_variable_slices[dimension_index].start,
                                              piece_read_slices[dimension_index].stop - input_variable_slices[dimension_index].start
                                              ) # No mask in this dimension
                                    )
                                    for dimension_index in range(len(input_variable.dimensions))
                                    ])
                                
                                logger.debug('piece_write_slices = {}'.format(piece_write_slices))
                                
                                # Get mask subsets for piece
                                piece_dim_masks = tuple([
                                    (variable_masks[dimension_index][piece_read_slices[dimension_index]] 
                                     if variable_masks[dimension_index] is not None else None)
                                    for dimension_index in range(len(input_variable.dimensions))
                                    ])
                                               
                                logger.debug('piece_dim_masks = {}'.format(pformat(piece_dim_masks)))
                                
                                # N.B: Nones in piece_dim_mask for unmasked dimensions will result in newaxis dimensions, but shape doesn't matter                               
                                output_variable[piece_write_slices] = input_variable[piece_read_slices][piece_dim_masks]
                                #logger.debug('output_variable[piece_write_slices] = {}'.format(output_variable[piece_write_slices]))
                        
                    else: # scalar variable - simple copy
                        logger.debug('\tCopying %s scalar data' % variable_name)
                        output_variable = input_variable
                else:
                    logger.debug('\tNot copying data for variable %s' % variable_name)
                    
            # Copy global attributes  
            logger.debug("Copying global attributes: %s" % ', '.join(self.netcdf_dataset.__dict__.keys()))
            for item, value in self.netcdf_dataset.__dict__.items():
                if type(value) == str:
                    nc_output_dataset.__setattr__(item, value.encode('utf-8'))
                else:
                    nc_output_dataset.__setattr__(item, value)
                    
            logger.debug('Finished copying netCDF dataset %s to %s.' % (self.nc_path, nc_out_path))
        
        finally:
            self.netcdf_dataset.set_auto_mask(True)
            nc_output_dataset.close()
            
    def get_crs_attributes(self, crs):
        '''\
        Function to return name and attributes of crs or transverse_mercator variable
        '''
        # Determine wkt and spatial_ref from crs as required
        if type(crs) == osgeo.osr.SpatialReference:
            spatial_ref = crs
        elif type(crs) == str:
            spatial_ref = get_spatial_ref_from_wkt(crs)
            
        wkt = spatial_ref.ExportToWkt() # Export WKT from spatial_ref for consistency even if supplied
        crs_attributes = {'spatial_ref': wkt}
        
        crs_attributes['inverse_flattening'] = spatial_ref.GetInvFlattening()
        crs_attributes['semi_major_axis'] = spatial_ref.GetSemiMajor()
        crs_attributes['longitude_of_prime_meridian'] = spatial_ref.GetAttrValue('PRIMEM', 1)

        if spatial_ref.GetUTMZone(): # CRS is UTM
            crs_variable_name = 'transverse_mercator'
            crs_attributes['grid_mapping_name'] = 'transverse_mercator'

            crs_attributes['latitude_of_projection_origin'] = spatial_ref.GetProjParm('latitude_of_origin')
            crs_attributes['scale_factor_at_central_meridian'] = spatial_ref.GetProjParm('scale_factor')
            crs_attributes['longitude_of_central_meridian'] = spatial_ref.GetProjParm('central_meridian')
            crs_attributes['false_northing'] = spatial_ref.GetProjParm('false_northing')
            crs_attributes['false_easting'] = spatial_ref.GetProjParm('false_easting')
                   
        else:
            #===================================================================
            # # Example crs attributes created by GDAL:
            #
            # crs:inverse_flattening = 298.257222101 ;
            # crs:spatial_ref = "GEOGCS[\"GEOCENTRIC DATUM of AUSTRALIA\",DATUM[\"GDA94\",SPHEROID[\"GRS80\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]" ;
            # crs:semi_major_axis = 6378137. ;
            # crs:GeoTransform = "121.1202390605822 0.0037 0 -20.56227098919639 0 -0.0037 " ;
            # crs:grid_mapping_name = "latitude_longitude" ;
            # crs:longitude_of_prime_meridian = 0. ;
            #===================================================================
            crs_variable_name = 'crs'
            crs_attributes['grid_mapping_name'] = 'latitude_longitude'

        logger.debug('{} attributes: {}'.format(pformat(crs_variable_name, crs_attributes)))
        return crs_variable_name, crs_attributes
    
    @abc.abstractmethod
    def get_convex_hull(self, to_wkt=None):
        '''\
        Abstract base function to return n x 2 array of coordinates for convex hull of all points
        Needs to be implemented in subclass (e.g. NetCDFPointUtils, NetCDFLineUtils, or NetCDFGridUtils)
        @param to_wkt: CRS WKT for shape
        '''
        pass
    
    @abc.abstractmethod
    def get_concave_hull(self, to_wkt=None, buffer_distance=None, tolerance=None):
        """\
        Abstract base function to return a shapely polygon for concave hull of all points
        Needs to be implemented in subclass (e.g. NetCDFPointUtils, NetCDFLineUtils, or NetCDFGridUtils)
        @param to_wkt: CRS WKT for shape
        @param buffer_distance: distance to buffer (kerf) initial shape outwards then inwards to simplify it
        @param tolerance: tolerance for simplification
        """
        pass
        
    def close(self):
        '''
        Function to close netCDF dataset if opened
        '''
        if self._netcdf_dataset:
            try:
                self._netcdf_dataset.close()
            except Exception as e:
                logger.warning('Unable to close {}: {}'.format(self.netcdf_path, e))
                
            self._netcdf_dataset = None
           
    @property
    def netcdf_dataset(self):
        '''
        Property getter function to open netCDF dataset only when required
        '''
        if not self._netcdf_dataset:
            logger.debug('Opening netCDF dataset {}'.format(self.nc_path))
            if self.opendap:
                try:
                    self._netcdf_dataset = netCDF4.Dataset(self.nc_path, mode="r")
                except OSError:
                    self._netcdf_dataset = netCDF4.Dataset(self.nc_path + '#fillmismatch', mode="r") # Work-around for _FillValue mismatch: https://github.com/Unidata/netcdf-c/issues/1299
            else:
                self._netcdf_dataset = netCDF4.Dataset(self.nc_path, mode="r")

        return self._netcdf_dataset
    
    @property
    def data_variable_list(self):
        '''
        Property getter function to return data_variable_list as required
        '''
        return self._data_variable_list

    @property
    def crs_variable(self):
        '''
        Property getter function to return crs_variable as required
        '''
        if self._crs_variable is None:
            logger.debug('Setting crs_variable property')
            for crs_variable_name in NetCDFUtils.CRS_VARIABLE_NAMES:
                self._crs_variable = self.netcdf_dataset.variables.get(crs_variable_name)
                
                if self._crs_variable is not None:
                    break
                
            #assert self._crs_variable is not None, 'Unable to determine crs_variable'
                
        return self._crs_variable


    @property
    def wkt(self):
        '''
        Property getter function to return wkt as required
        '''
        if not self._wkt:
            logger.debug('Setting wkt property')
            self._wkt = self.crs_variable.spatial_ref
        return self._wkt


    @property
    def x_variable(self):
        return self._x_variable
    
    @property
    def y_variable(self):
        return self._y_variable
    
    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, debug_value):
        if self._debug != debug_value or self._debug is None:
            self._debug = debug_value
            if self._debug:
                log_level = logging.DEBUG
            else:
                log_level = logging.INFO
                
            logger.setLevel(log_level)
            logger.debug('Logger {} set to level {}'.format(logger.name, log_level))
                
            logging.getLogger(self.__module__).setLevel(log_level)
            logger.debug('Logger {} set to level {}'.format(self.__module__, log_level))
            
            for base_class in self.__class__.__bases__:
                logging.getLogger(base_class.__module__).setLevel(log_level)
                logger.debug('Logger {} set to level {}'.format(base_class.__module__, log_level))

                
    @property
    def wgs84_bbox(self):
        '''
        Property getter function to return wgs84_bbox as required
        '''
        if self._wgs84_bbox is None:
            logger.debug('Setting wgs84_bbox property')
            transformed_bbox = transform_coords(self.native_bbox, from_wkt=self.wkt, to_wkt='EPSG:4326')
            self._wgs84_bbox = [[min(transformed_bbox[:,0]), min(transformed_bbox[:,1])],
                                [max(transformed_bbox[:,0]), min(transformed_bbox[:,1])],
                                [max(transformed_bbox[:,0]), max(transformed_bbox[:,1])],
                                [min(transformed_bbox[:,0]), max(transformed_bbox[:,1])]
                                ]
        return self._wgs84_bbox


def main():
    '''
    Main function for calling NetCDFUtils.copy function
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-f", "--format", help="NetCDF file format (one of 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA')",
                        type=str, default='NETCDF4')
    parser.add_argument("--chunkspec", help="comma-separated list of <dimension_name>/<chunk_size> specifications",
                        type=str)
    parser.add_argument("--complevel", help="Compression level for chunked variables as an integer 0-9. Default is 4",
                        type=int, default=4)
    parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    
    args = parser.parse_args()
    
    if args.chunkspec:
        chunk_spec = {dim_name: int(chunk_size) 
                    for dim_name, chunk_size in [chunk_spec_string.strip().split('/') for chunk_spec_string in args.chunkspec.split(',')]}
    else:
        chunk_spec = None
            
    ncu = NetCDFUtils(args.input_path,
                      debug=args.debug
                      )   
    ncu.max_bytes = 50000
    print('max_bytes = {}'.format(ncu.max_bytes))
          
    ncu.copy(args.output_path, 
             #datatype_map_dict={},
             # Compress all chunked variables
             variable_options_dict={variable_name: {'chunksizes': [chunk_spec.get(dimension) 
                                                                   for dimension in variable.dimensions
                                                                   ],
                                                    'zlib': bool(args.complevel),
                                                    'complevel': args.complevel
                                                    }
                               for variable_name, variable in ncu.netcdf_dataset.variables.items()
                               if (set(variable.dimensions) & set(chunk_spec.keys()))
                               } if chunk_spec else {},
             #dim_range_dict={'lat': (5,205),'lon': (5,305)},
             #dim_mask_dict={},
             nc_format=args.format,
             #limit_dim_size=False
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
        

