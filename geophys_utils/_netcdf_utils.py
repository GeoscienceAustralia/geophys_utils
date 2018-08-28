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

import netCDF4
import math
import itertools
import argparse
import re
from distutils.util import strtobool
from geophys_utils._crs_utils import get_spatial_ref_from_wkt
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

class NetCDFUtils(object):
    '''
    NetCDFUtils class implementing useful functionality against netCDF files
    '''
    DEFAULT_COPY_OPTIONS = {'complevel': 4, 
                            'zlib': True, 
                            'fletcher32': True,
                            'shuffle': True,
                            'endian': 'little',
                            'chunksizes': [1024, 1024]
                            }

    def __init__(self, netcdf_dataset, debug=False):
        '''
        Constructor for NetCDFUtils
        '''
        self._debug = None # Initialise private variable
        self.debug = debug # Set debug property
        
        if type(netcdf_dataset) == str: # String provided as path to netCDF file
            self.nc_path = netcdf_dataset
            self.netcdf_dataset = netCDF4.Dataset(self.nc_path, mode="r")
        elif type(netcdf_dataset) == netCDF4.Dataset: # NetCDF4.Dataset object provided
            self.netcdf_dataset = netcdf_dataset 
            self.nc_path = netcdf_dataset.filepath()
        
        self.opendap = (re.match('^http.*', self.nc_path) is not None)
        if self.opendap:
            self.max_bytes = 500000000 # 500MB limit for NCI OPeNDAP
        else:
            self.max_bytes = 4000000000 # 4GB limit for direct netCDF file access
        
        # Identify all spatial grid variables by first finding the variable with 'grid_mapping'
        # and then finding all variables with the same dimensionality
        #TODO: Do something not specific to gridded datasets
        self.data_variable_list = []
        try:
            data_variable_dimensions = [variable for variable in self.netcdf_dataset.variables.values() 
                                       if hasattr(variable, 'grid_mapping')][0].dimensions
            self.data_variable_list = [variable for variable in self.netcdf_dataset.variables.values() 
                                       if variable.dimensions == data_variable_dimensions]
        except:
            logger.debug('Unable to determine data variable(s) (must have same dimensions as variable with "grid_mapping" attribute)')

        
        #TODO: Make sure this is general for all CRSs
        self.x_variable = (self.netcdf_dataset.variables.get('lon') 
                           or self.netcdf_dataset.variables.get('x')
                           )
        
        self.y_variable = (self.netcdf_dataset.variables.get('lat') 
                           or self.netcdf_dataset.variables.get('y')
                           )
        
        self.y_inverted = (self.y_variable[-1] < self.y_variable[0]) if self.y_variable else False
        
        try:
            self.grid_mapping_variable = self.netcdf_dataset.variables[
                self.data_variable_list[0].grid_mapping]
        except:
            self.grid_mapping_variable = None
            for grid_mapping_variable_name in ['crs',
                                               'transverse_mercator'
                                               ]:
                try:
                    self.grid_mapping_variable = self.netcdf_dataset.variables[grid_mapping_variable_name]
                    break
                except: 
                    continue

            
        try:
            self.wkt = self.grid_mapping_variable.spatial_ref
        except:
            try:
                self.wkt = get_spatial_ref_from_wkt(self.grid_mapping_variable.epsg_code).ExportToWkt()
            except:
                #TODO: Do something a bit better than assuming unprojected WGS84
                self.wkt = get_spatial_ref_from_wkt('EPSG:4326').ExportToWkt()

    def copy(self, nc_out_path, 
                 datatype_map_dict={},
                 variable_options_dict={},
                 dim_range_dict={},
                 nc_format=None,
                 limit_dim_size=False,
                 empty_var_list=[],
                 invert_y=None):
        '''
        Function to copy a netCDF dataset to another one with potential changes to size, format, 
            variable creation options and datatypes.
            
            @param nc_in_path: path to existing netCDF input file 
            @param nc_out_path: path to netCDF output file 
            @param datatype_map_dict: dict containing any maps from source datatype to new datatype.
                e.g. datatype_map_dict={'uint64': 'uint32'}  would convert all uint64 variables to uint32.
            @param variable_options_dict: dict containing any overrides for per-variable variable creation 
                options. e.g. variable_options_dict={'sst': {'complevel': 2, 'zlib': True}} would apply
                compression to variable 'sst'
            @param dim_range_dict: dict of (start, end+1) tuples keyed by dimension name
            @param nc_format: output netCDF format - 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET', 
                'NETCDF3_64BIT_DATA', 'NETCDF4_CLASSIC', or 'NETCDF4'. Defaults to same as input format.  
            @param limit_dim_size: Boolean flag indicating whether unlimited dimensions should be fixed
            @param empty_var_list: List of strings denoting variable names for variables which should be created but not copied
            @param invert_y: Boolean parameter indicating whether copied Y axis should be Southwards positive (None means same as source)
        '''  
             
        if invert_y is None: # Default y-axis inversion to same as source
            invert_y = self.y_inverted
            
        flip_y = (invert_y != self.y_inverted)
        
        # Override default variable options with supplied ones for all data variables
        for data_variable in self.data_variable_list:
            variable_dict = dict(NetCDFUtils.DEFAULT_COPY_OPTIONS)
            variable_dict.update(variable_options_dict.get(data_variable.name) or {})
            variable_options_dict[data_variable.name] = variable_dict
                                
        nc_format = nc_format or self.netcdf_dataset.file_format 
        logger.info('Output format is %s' % nc_format)
        
        nc_output_dataset = netCDF4.Dataset(nc_out_path, mode="w", clobber=True, format=nc_format)
        
        try:
            dims_used = set()
            dim_size = {}
            for variable_name, variable in self.netcdf_dataset.variables.items():
                dims_used |= set(variable.dimensions)
                
                for dimension_index in range(len(variable.dimensions)): 
                    dimension_name = variable.dimensions[dimension_index]
                    if dim_size.get(dimension_name):
                        continue
                        
                    dim_range = dim_range_dict.get(dimension_name)
                    if dim_range:
                        dim_size[dimension_name] = dim_range[1] - dim_range[0]
                    else:
                        dim_size[dimension_name] = variable.shape[dimension_index]
                    
            #logger.debug(dim_size)
            
            #Copy dimensions
            for dimension_name, dimension in self.netcdf_dataset.dimensions.items():
                if dimension_name in dims_used: # Discard unused dimensions
                    logger.info('Copying dimension %s of length %d' % (dimension_name, dim_size[dimension_name]))
                    nc_output_dataset.createDimension(dimension_name, 
                                          dim_size[dimension_name] 
                                          if not dimension.isunlimited() or limit_dim_size 
                                          else None)
                else:
                    logger.info('Skipping unused dimension %s' % dimension_name)
    
            # Copy variables
            for variable_name, input_variable in self.netcdf_dataset.variables.items():
                dtype = datatype_map_dict.get(str(input_variable.datatype)) or input_variable.datatype
                
                # Special case for "crs" or "transverse_mercator" - want byte datatype
                if input_variable == self.grid_mapping_variable: 
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
                        var_options['chunksizes'][dimension_index] = min(var_options['chunksizes'][dimension_index],
                                                                         dim_size[input_variable.dimensions[dimension_index]])
                             
                options_string = ' with options: %s' % ', '.join(['%s=%s' % item for item in var_options.items()]) if var_options else ''   
                logger.info("Copying variable %s from datatype %s to datatype %s%s" % (variable_name, 
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
                logger.info('\tCopying %s attributes: %s' % (variable_name, ', '.join(input_variable.ncattrs())))
                output_variable.setncatts({k: input_variable.getncattr(k) for k in input_variable.ncattrs() if not k.startswith('_')})
                
                #===============================================================
                # if (flip_y and (input_variable == self.grid_mapping_variable)):                    
                #     output_GeoTransform = list(self.GeoTransform)
                #     output_GeoTransform[5] = - output_GeoTransform[5]
                #     output_variable.GeoTransform = ' '.join([str(value) for value in output_GeoTransform])
                #     logger.info('%s.GeoTransform rewritten as "%s"' % (variable_name, output_variable.GeoTransform))
                #===============================================================
    
                if variable_name not in empty_var_list:
                    # Copy data
                    if input_variable.shape: # array
                        overall_slices = [slice(*dim_range_dict[input_variable.dimensions[dimension_index]])  
                                  if dim_range_dict.get(input_variable.dimensions[dimension_index])
                                  else slice(0, input_variable.shape[dimension_index])
                                  for dimension_index in range(len(input_variable.dimensions))
                                 ]
                        #logger.debug('overall_slices={}.format(overall_slices))
                        logger.info('\tCopying %s array data of shape %s' % (variable_name,
                                                                             tuple([overall_slices[dimension_index].stop - overall_slices[dimension_index].start
                                                                                    for dimension_index in range(len(input_variable.dimensions))]
                                                                                   )
                                                                             )
                                    )
                        
                        if (not input_variable_chunking or 
                            len(input_variable.dimensions) != 2): 
                            # No chunking - Try to copy in one hit
                            
                            if ((input_variable == self.y_variable) and flip_y): 
                                # Y-axis flip required
                                assert len(overall_slices) == 1, 'y-axis variable should be one-dimensional'
                                overall_slices = [slice(overall_slices[0].stop-1, overall_slices[0].start-1 if overall_slices[0].start else None, -1)]
                                logger.info('\tInverting y-axis variable %s' % variable_name)
                                
                            output_variable[...] = input_variable[overall_slices]
                        
                        else: # Chunked - perform copy in pieces
                            #TODO: Improve this for small chunks
                            assert len(input_variable.dimensions) == 2, 'Can only chunk copy 2D data at the moment'
                            
                            # Use largest chunk sizes between input and output
                            piece_sizes = [max(var_options['chunksizes'][dimension_index],
                                              input_variable_chunking[dimension_index])
                                          for dimension_index in range(len(input_variable.dimensions))
                                         ]
        
                            piece_index_ranges = [(overall_slices[dimension_index].start // piece_sizes[dimension_index],
                                                  int(math.ceil(float(overall_slices[dimension_index].stop) / piece_sizes[dimension_index]))
                                                 )
                                                 for dimension_index in range(len(input_variable.dimensions))
                                                ]
                                                 
                            piece_counts = [int(math.ceil(float(dim_size[input_variable.dimensions[dimension_index]]) / 
                                                          piece_sizes[dimension_index]))
                                            for dimension_index in range(len(input_variable.dimensions)) 
                                            ]
                        
                            logger.info('\tCopying %s pieces of size %s cells' % (' x '.join([str(piece_count) for piece_count in piece_counts]),
                                                                          ' x '.join([str(piece_size) for piece_size in piece_sizes])
                                                                          )
                                        )
                            
                            try:
                                ydim_index = input_variable.dimensions.index(self.y_variable.name)
                            except ValueError:
                                ydim_index = None
                            
                            # Iterate over every piece
                            for piece_indices in itertools.product(*[range(piece_index_ranges[dimension_index][0], 
                                                                          piece_index_ranges[dimension_index][1])
                                                               for dimension_index in range(len(input_variable.dimensions))
                                                              ]
                                                            ):
                                
                                                   
                                logger.info('\t\tCopying piece %s' % (piece_indices,))
                                
                                piece_read_slices = [slice(max(overall_slices[dimension_index].start,
                                                         piece_indices[dimension_index] * piece_sizes[dimension_index]
                                                        ),                                                 
                                                     min(overall_slices[dimension_index].stop,
                                                         (piece_indices[dimension_index] + 1) * piece_sizes[dimension_index]
                                                        )
                                                    )
                                               for dimension_index in range(len(input_variable.dimensions))
                                               ]
                                
                                piece_write_slices = [slice(piece_read_slices[dimension_index].start - overall_slices[dimension_index].start,
                                                      piece_read_slices[dimension_index].stop - overall_slices[dimension_index].start,
                                                     )
                                                for dimension_index in range(len(input_variable.dimensions))
                                               ]
                                
                                if flip_y and ydim_index is not None:
                                    # Flip required
                                    piece_write_slices[ydim_index] = slice(output_variable.shape[ydim_index] - piece_write_slices[ydim_index].start -1, 
                                                                          output_variable.shape[ydim_index] - piece_write_slices[ydim_index].stop - 1 
                                                                            if (output_variable.shape[ydim_index] - piece_write_slices[ydim_index].stop) 
                                                                            else None, -1)
                                
                                #logger.debug(piece_read_slices, piece_write_slices)
                                
                                output_variable[piece_write_slices] = input_variable[piece_read_slices]
                        
                    else: # scalar variable - simple copy
                        logger.info('\tCopying %s scalar data' % variable_name)
                        output_variable = input_variable
                else:
                    logger.info('\tNot copying data for variable %s' % variable_name)
                    
            # Copy global attributes  
            logger.info("Copying global attributes: %s" % ', '.join(self.netcdf_dataset.__dict__.keys()))
            for item, value in self.netcdf_dataset.__dict__.items():
                if type(value) == str:
                    nc_output_dataset.__setattr__(item, value.encode('utf-8'))
                else:
                    nc_output_dataset.__setattr__(item, value)
                    
            logger.info('Finished copying netCDF dataset %s to %s.' % (self.nc_path, nc_out_path))
        
        finally:
            nc_output_dataset.close()
            
            
    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, debug_value):
        if self._debug != debug_value or self._debug is None:
            self._debug = debug_value
            
            if self._debug:
                logger.setLevel(logging.DEBUG)
                logging.getLogger(self.__module__).setLevel(logging.DEBUG)
            else:
                logger.setLevel(logging.INFO)
                logging.getLogger(self.__module__).setLevel(logging.INFO)
                
        logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))
        logging.getLogger(self.__module__).debug('Logger {} set to level {}'.format(self.__module__, logger.level))

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
    parser.add_argument("--chunking", help="comma-separated list of chunk sizes for each dimension",
                        type=str)
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    parser.add_argument('-i', '--invert_y', help='Store copy with y-axis indexing Southward positive', type=str)
    
    args = parser.parse_args()
    
    if args.invert_y is not None:
        invert_y = bool(strtobool(args.invert_y))
    else:
        invert_y = None # Default to same as source
    
    if args.do_copy:
        if args.chunking:
            chunking = [int(chunk_size.strip()) for chunk_size in args.chunking.split(',')]
        else:
            chunking = None
            
    ncu = NetCDFUtils(args.input_path)   
    
    ncu.copy(args.output_path, 
             #datatype_map_dict={},
             variable_options_dict={data_variable.name: {'chunksizes': chunking}
                               for data_variable in ncu.data_variable_list
                               } if chunking else {},
             #dim_range_dict={},
             #nc_format=None,
             #limit_dim_size=False
             invert_y=invert_y
             )
        

if __name__ == '__main__':
    main()
