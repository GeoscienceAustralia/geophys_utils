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
ASEGGDF2NetCDFConverter concrete class for converting ASEG-GDF data to netCDF

Created on 28Mar.2018

@author: Alex Ip
'''
import argparse
from collections import OrderedDict
import numpy as np
import re
import os
import sys
from datetime import datetime
from pprint import pformat
import yaml
import tempfile
import netCDF4
import logging

from geophys_utils.netcdf_converter import ToNetCDFConverter, NetCDFVariable
from geophys_utils import get_spatial_ref_from_wkt
from geophys_utils.netcdf_converter.aseg_gdf_utils import aseg_gdf_format2dtype, fix_field_precision, truncate
from geophys_utils import points2convex_hull
from geophys_utils import transform_coords

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

TEMP_DIR = tempfile.gettempdir()
#TEMP_DIR = 'D:\Temp'
#TEMP_DIR = 'U:\Alex\Temp'

# Set this to zero for no limit - only set a non-zero value for testing
POINT_LIMIT = 0

# Number of rows per chunk in temporary netCDF cache file
CACHE_CHUNK_ROWS = 8192

class ASEGGDF2NetCDFConverter(ToNetCDFConverter):
    '''
    ASEGGDF2NetCDFConverter concrete class for converting ASEG-GDF data to netCDF
    '''
    def __init__(self, 
                 nc_out_path, 
                 aem_dat_path, 
                 dfn_path, 
                 crs_string=None,
                 netcdf_format='NETCDF4', 
                 default_chunk_size=None, 
                 default_variable_parameters=None,
                 settings_path=None,
                 fix_precision=True,
                 space_delimited=False
                 ):
        '''
        Concrete constructor for subclass ASEGGDF2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure the base class constructor is called from the subclass constructor
        @param nc_out_path: Path to output netCDF file on filesystem
        @param aem_dat_path: Path to .dat AEM data source file on filesystem
        @param dfn_path: Path to .dfn definition file on filesystem
        @param netcdf_format: Format for netCDF file. Defaults to 'NETCDF4_CLASSIC'
        @param default_chunk_size: single default chunk size for all dimensions. None means take default, zero means not chunked.
        @param default_variable_parameters: Optional dict containing default parameters for netCDF variable creation
        @param settings_path: Optional path for settings YAML file
        @param fix_precision: Optional Boolean flag indicating whether to fix (i.e. reduce) field precisions
        '''

        def get_field_definitions():
            '''
            Function to read raw field definitions from .dfn file
            Will set self.dimensions as an Ordereddict of dimension sizes keyed by dimension name
            '''
            def parse_dfn_file(dfn_path):
                logger.info('Reading definitions file {}'.format(dfn_path))
                self.field_definitions = []
                
                dfn_file = open(dfn_path, 'r')
                for line in dfn_file:
                    key_value_pairs = {}
                    positional_value_list = []
                    for semicolon_split_string in [semicolon_split_string.strip() for semicolon_split_string in line.split(';')]:
                        for colon_split_string in [colon_split_string.strip() for colon_split_string in semicolon_split_string.split(':')]:
                            for comma_split_string in [comma_split_string.strip() for comma_split_string in colon_split_string.split(',')]:
                                definition = [equals_split_string.strip() for equals_split_string in comma_split_string.split('=')]
                                if len(definition) == 2:
                                    key_value_pairs[definition[0]] = definition[1]
                                elif len(definition) == 1:
                                    positional_value_list.append(definition[0])
    
                    logger.debug('key_value_pairs: {},\npositional_value_list: {}'.format(pformat(key_value_pairs), pformat(positional_value_list))) 
                
                    # Column definition
                    if key_value_pairs.get('RT') in ['', 'DATA'] and (positional_value_list 
                                                            and positional_value_list[0] != 'END DEFN'): 
                        short_name = positional_value_list[0].lower()
                        fmt = positional_value_list[1] if len(positional_value_list) >= 2 else None
                        units = key_value_pairs.get('UNITS') or key_value_pairs.get('UNIT')
                        long_name = key_value_pairs.get('NAME') or (positional_value_list[2] if len(positional_value_list) >= 3 else None)
                        fill_value = float(key_value_pairs.get('NULL')) if key_value_pairs.get('NULL') is not None else None
                        
                        # Parse format to determine columns, data type and numeric format
                        dtype, columns, width_specifier, decimal_places = aseg_gdf_format2dtype(fmt)
                                    
                        field_definition = {'short_name': short_name,
                                      'format': fmt,
                                      'long_name': long_name,
                                      'dtype': dtype,
                                      'columns': columns,
                                      'width_specifier': width_specifier, 
                                      'decimal_places': decimal_places
                                      }
                        if units:
                            field_definition['units'] = units
                        if fill_value is not None:
                            field_definition['fill_value'] = fill_value
                             
                        
                        # Set variable attributes in field definition
                        variable_attribute_dict = {attribute_name: key_value_pairs.get(key.upper())
                            for key, attribute_name in self.settings['variable_attributes'].items()
                            if key_value_pairs.get(key.upper()) is not None
                                                   }
                        
                        # Store aseg_gdf_format in variable attributes
                        variable_attribute_dict['aseg_gdf_format'] = fmt 

                        if variable_attribute_dict:
                            field_definition['variable_attributes'] = variable_attribute_dict    
                        
                        self.field_definitions.append(field_definition)
                        
                                    
                    # Set CRS from projection name
                    elif not self.spatial_ref:
                        if (key_value_pairs.get('RT') == 'PROJ') and (positional_value_list[0] == 'COORDSYS'): # As per ASEG standard                       
                            projection_name = key_value_pairs.get('NAME')
                            if projection_name:
                                self.spatial_ref = get_spatial_ref_from_wkt(projection_name)
                                logger.debug('CRS set from .dfn file COORDSYS COMMENT attribute {}'.format(projection_name))
                                break # Nothing more to do
                        elif (key_value_pairs.get('RT') == 'PROJ') and (positional_value_list[0] == 'PROJNAME'): # Non-standard                        
                            projection_name = key_value_pairs.get('COMMENT')
                            if projection_name:
                                self.spatial_ref = get_spatial_ref_from_wkt(projection_name)
                                logger.debug('CRS set from .dfn file PROJNAME NAME attribute {}'.format(projection_name))
                                break # Nothing more to do
                        elif (key_value_pairs.get('RT') == 'PROJ') and (positional_value_list[0] == 'DATUM'): # Unprojected                         
                            projection_name = key_value_pairs.get('NAME')
                            if projection_name:
                                self.spatial_ref = get_spatial_ref_from_wkt(projection_name)
                                logger.debug('CRS set from .dfn file DATUM NAME attribute {}'.format(projection_name))
                                break # Nothing more to do
                            
            # Start of get_field_definitions function
            parse_dfn_file(dfn_path)
            
            # Read overriding field definition values from settings
            if self.settings.get('field_definitions'):
                for field_definition in self.field_definitions:
                    overriding_field_definition = self.settings['field_definitions'].get(field_definition['short_name'])
                    if overriding_field_definition:
                        field_definition.update(overriding_field_definition)
            
            logger.debug('self.dimensions: {}'.format(pformat(self.dimensions)))
            logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
                        
            # Check for CRS definition in latitude field - need to do this after short_name has potentially been remapped
            if not self.spatial_ref:
                try:
                    field_definition = [field_definition for field_definition in self.field_definitions if field_definition['short_name'] == 'latitude'][0]
                    crs_string = field_definition['variable_attributes']['datum_name']
                    self.spatial_ref = get_spatial_ref_from_wkt(crs_string)
                    logger.debug('CRS set from latitude variable datum_name attribute {}'.format(crs_string))
                except:
                    logger.debug('Unable to set CRS from latitude datum_name attribute')               

            assert self.spatial_ref, 'Coordinate Reference System undefined'
            logger.debug('self.spatial_ref: {}'.format(self.spatial_ref.ExportToPrettyWkt()))
        
        def read_data_file(): 
            '''
            Function to read data file into temporary netCDF cache
            '''   
            def create_nc_cache():
                '''
                Function to create temporary cache file with one 2D variable
                Needs to have self.column_count defined to work
                '''
                self.nc_cache_path = os.path.join(TEMP_DIR, re.sub('\W+', '_', os.path.splitext(self.aem_dat_path)[0]) + '.nc')
                self._nc_cache_dataset = netCDF4.Dataset(self.nc_cache_path, mode="w", clobber=True, format='NETCDF4')
                self._nc_cache_dataset.createDimension(dimname='rows', size=None) # Unlimited size
                
                self.column_count = 0
                for field_definition in self.field_definitions:
                    short_name = field_definition['short_name']
                    columns = field_definition['columns']

                    if columns == 1: # 1D variable
                        field_dimensions =('rows',)
                        chunksizes=(CACHE_CHUNK_ROWS,)
                    else: # 2D variable
                        field_dimension_name = short_name + '_dim' # Default 2nd dimension name

                        # Look up any dimension name(s) for this variable from settings
                        override_definition = self.settings['field_definitions'].get(short_name)
                        # Determine name of dimension
                        if override_definition:
                            override_dimensions = override_definition.get('dimensions')
                            if override_dimensions:
                                if type(override_dimensions) == list:
                                    field_dimension_name = '_plus_'.join(override_dimensions) # Multiple partial dimensions
                                elif type(override_dimensions) == str:
                                    # Note that the "_dim" suffix is required to avoid triggering an obscure netCDF4 bug which seems to occur when 
                                    # a dimension is defined with the same name as an existing variable.
                                    # See https://github.com/Unidata/netcdf-c/issues/295
                                    # This does not appear to be an issue when the dimension is defined before the variable, but we can't do that for
                                    # the cache dataset, only the final one
                                    field_dimension_name = override_dimensions + '_dim' # Single dimension
                                    
                        # Create dimension if it doesn't already exist
                        if field_dimension_name not in self._nc_cache_dataset.dimensions.keys():
                            self._nc_cache_dataset.createDimension(dimname=field_dimension_name, size=columns)
                            
                        field_dimensions =('rows', field_dimension_name)
                        chunksizes=(CACHE_CHUNK_ROWS, columns)
                        
                    logger.debug('\tCreating cache variable {} with dtype {} and dimension(s) {}'.format(short_name, field_definition['dtype'], field_dimensions))
                    self._nc_cache_dataset.createVariable(varname=short_name,
                                                          datatype=field_definition['dtype'],
                                                          dimensions=field_dimensions,
                                                          chunksizes=chunksizes,
                                                          **NetCDFVariable.DEFAULT_VARIABLE_PARAMETERS
                                                          )
                    self.column_count += columns
                    
                try: # Do a sync now to force an error if there are dimension/variable naming conflicts
                    self._nc_cache_dataset.sync()
                    logger.debug('NetCDF sync completed')
                except:
                    logger.debug(self._nc_cache_dataset.dimensions, self._nc_cache_dataset.variables)
                    raise
                              
                logger.debug('Created temporary cache file {}'.format(self.nc_cache_path))
            
            def cache_chunk_list(chunk_list, start_row):
                '''
                Helper function to write list of lists to cache variables
                '''
                end_row = start_row + len(chunk_list)
                logger.debug('Writing rows {}-{} to disk cache'.format(start_row, end_row))
                for field_index in range(len(self.field_definitions)):
                    field_definition = self.field_definitions[field_index]
                    short_name = field_definition['short_name']
                    cache_variable = self._nc_cache_dataset.variables[short_name]
                    
                    chunk_array = np.array([row_list[field_index] for row_list in chunk_list])
                    #logger.debug('{} cache_variable: {}'.format(cache_variable.name, cache_variable))
                    #logger.debug('{}: {}-element {} chunk_array: {}'.format(short_name, chunk_array.shape[0], chunk_array.dtype, chunk_array))

                    cache_variable[start_row:end_row] = chunk_array
            
                self.total_points += len(chunk_list)
                
                
            def read_fixed_length_fields(line):
                '''
                Helper function to read fixed length fields into a list of lists
                '''
                row_list = []
                line_column_count = 0
                start_char = 0
                for field_definition in self.field_definitions:
                    short_name = field_definition['short_name']
                    columns = field_definition['columns']
                    dtype = field_definition['dtype']
                    aseg_gdf_format = field_definition['format']
                    column_list = []
                    for _column_offset in range(columns):
                        end_char = start_char + field_definition['width_specifier']
                        value_string = line[start_char:end_char]
                        
                        # Work-around for badly formatted files with first entry too short
                        if not aseg_gdf_format.startswith('A') and ' ' in value_string.strip(): # Not a string field and has a space in the middle
                            value_string = re.match('\s*\S*', value_string).group(0) # Strip anything after non-leading whitespace character
                            end_char = start_char + len(value_string) # Adjust character offset for next column
                        
                        value_string = value_string.strip()
                        try:
                            if dtype.startswith('int'):
                                value = int(value_string)
                            elif dtype.startswith('float'):
                                value = float(value_string)
                            else: # Assume string
                                value = value_string
                        except ValueError:
                            logger.warning('Unable to convert "{}" field value "{}" to type {}'.format(short_name, value_string, dtype))
                            logger.debug('line: "{}"'.format(line))
                            return None
                        #logger.debug('value_string: {}, repr(value): {}'.format(value_string, repr(value)))
                        line_column_count += 1
                        column_list.append(value)
                        start_char = end_char
                    #logger.debug('column_list: {}'.format(column_list))
                    if column_list:
                        if columns == 1: # 1D variable
                            row_list.append(column_list[0]) 
                        else:
                            row_list.append(column_list) 
                            
                #logger.debug('row_list: {}'.format(row_list))
                if line_column_count != self.column_count:
                    logger.warning('Invalid number of columns found: Expected {}, found {}. Skipping line'.format(self.column_count, line_column_count))
                    logger.debug('line: "{}"'.format(line))
                    return None
                return row_list
                
                
            def read_delimited_fields(line, delimiter_regex='\s+'):
                '''
                Helper function to read delimited fields into a list of lists
                N.B: This should not be required, but ASEG-supplied example file Example_Gravity_Springfield_1989.dat requires it
                '''
                row_list = []
                value_list = re.sub(delimiter_regex, '\t', line.strip()).split('\t')

                #logger.debug('value_list: {}'.format(value_list))
                #logger.debug('row_list: {}'.format(row_list))
                if len(value_list) != self.column_count:
                    logger.warning('Invalid number of columns found: Expected {}, found {}. Skipping line'.format(self.column_count, len(value_list)))
                    logger.debug('line: "{}"'.format(line))
                    return None

                value_index = 0
                for field_definition in self.field_definitions:
                    short_name = field_definition['short_name']
                    columns = field_definition['columns']
                    dtype = field_definition['dtype']
                    column_list = []
                    for _column_offset in range(columns):
                        value_string = value_list[value_index].strip()
                        try:
                            if dtype.startswith('int'):
                                value = int(value_string)
                            elif dtype.startswith('float'):
                                value = float(value_string)
                            else: # Assume string
                                value = value_string
                        except ValueError:
                            logger.warning('Unable to convert "{}" field value "{}" to type {}'.format(short_name, value_string, dtype))
                            logger.debug('line: "{}"'.format(line))
                            return None
                        #logger.debug('value_string: {}, repr(value): {}'.format(value_string, repr(value)))
                        value_index += 1
                        column_list.append(value)

                    #logger.debug('column_list: {}'.format(column_list))
                    if column_list:
                        if columns == 1: # 1D variable
                            row_list.append(column_list[0]) 
                        else:
                            row_list.append(column_list) 
                            
                return row_list
                
                
            logger.info('Reading data file {}'.format(aem_dat_path))
            aem_dat_file = open(aem_dat_path, 'r')
            
            create_nc_cache()
            
            line_count = 0
            self.total_points = 0
            chunk_list = []
            for line in aem_dat_file:
                line = re.sub('\n$', '', line) # Strip trailing EOL
                #logger.debug('line: "{}"'.format(line))
                if not line.strip(): # Skip empty lines
                    continue
                
                if self.space_delimited:
                    row_list = read_delimited_fields(line, '\s+')
                elif '\t' in line: # Assume tab delimited line
                    row_list = read_delimited_fields(line, '\t')
                else:
                    row_list = read_fixed_length_fields(line)
                    
                if not row_list: # Invalid line - ignore
                    continue
                
                chunk_list.append(row_list)
                line_count += 1
                
                if len(chunk_list) % CACHE_CHUNK_ROWS == 0:
                    cache_chunk_list(chunk_list, self.total_points)
                    chunk_list = [] # Reset chunk after writing
                
                if not line_count % 10000:
                    logger.info('{} lines read'.format(line_count))
                     
                if POINT_LIMIT and self.total_points >= POINT_LIMIT:
                    logger.debug('Truncating input for testing after {} points'.format(POINT_LIMIT))
                    break
                 
            if len(chunk_list):
                cache_chunk_list(chunk_list, self.total_points)
                chunk_list = [] # Reset chunk after writing
            
            aem_dat_file.close()             
    
            logger.info('A total of {} points were read'.format(self.total_points))
            
            assert self.total_points, 'Unable to read any points'
            
            self.dimensions['point'] = self.total_points # All files will have a point dimension - declare this first
            
            
        def modify_field_definitions():    
            '''
            Function to update field definitions and dimensions based on values read from data file
            Assumes all dimension fields declared before 2D fields
            '''
            column_start_index = 0
            field_definition_index = -1
            while field_definition_index < len(self.field_definitions) - 1:
                field_definition_index += 1
                field_definition = self.field_definitions[field_definition_index]
                field_name = field_definition['short_name']
                cache_variable = self._nc_cache_dataset.variables[field_name]
                logger.debug('field_definition_index: {}, field_name: {}'.format(field_definition_index, field_name))
                
                if field_name in self.settings['dimension_fields']:
                    self.dimensions[field_name] = int(cache_variable[0]) # Read value from first line
                    #logger.debug('self.dimensions["{}"] set to {}'.format(field_name, cache_variable[0]))
                    field_definition['column_start_index'] = None # Do not output this column
                    column_start_index += 1 # Single column only
                    continue # Check next field
                
    
                # Non-dimension field                
                field_definition['column_start_index'] = column_start_index
                
                # Check for 2D column count in format and strip integers from start of format if found
                match = re.match('^\d+', field_definition['format'])
                if match: # format string starts with integers - 2D variable
                    field_dimension_size = int(match.group(0))
                    field_definition['format'] = re.sub('^\d+', '', field_definition['format'])
                    default_dimension_names = field_definition.get('dimensions') # Look for default dimension name
                    #logger.debug('default_dimension_names: {}'.format(default_dimension_names))
                    #logger.debug('self.dimensions: {}'.format(self.dimensions))
                    if default_dimension_names: # Default dimension name(s) found from settings
                        if type(default_dimension_names) == str: # Only one dimension specified as a string
                            actual_dimension_size = self.dimensions.get(default_dimension_names)
                            logger.debug('actual_dimension_size: {}'.format(actual_dimension_size))
                            if actual_dimension_size is None: # Default dimension not already defined - create it
                                self.dimensions[default_dimension_names] = field_dimension_size
                                #logger.debug('Set self.dimensions[{}] = {}'.format(default_dimension_names, field_dimension_size))
                            else:
                                assert field_dimension_size == actual_dimension_size, 'Mismatched dimension sizes'
                                
                            logger.debug('Variable {} has dimension {} of size {}'.format(field_name,
                                                                                      default_dimension_names, 
                                                                                      field_dimension_size))
                        
                        # Use dimension size values read from data file to create multiple part-arrays from a single composite
                        # N.B: This could be considered to be encouraging bad behaviour - the .dfn file should 
                        # ideally be re-written
                        elif type(default_dimension_names) == list: # Multiple part-dimensions specified
                            assert set(default_dimension_names) <= set(self.dimensions.keys()), 'Dimensions {} not all defined in {}'.format(default_dimension_names, self.dimensions.keys())
                            combined_dimension_size = sum([self.dimensions[key] for key in default_dimension_names])
                            assert combined_dimension_size == field_dimension_size, 'Incorrect dimension sizes. Expected {}, found {}'.format(combined_dimension_size, field_dimension_size)
                            logger.debug('Removing original composite field definition')
                            self.field_definitions.remove(field_definition)
                            field_definition_index -= 1
                            # Create new field definitions for part dimensions
                            column_offset = 0
                            for part_dimension_name in default_dimension_names:
                                columns = self.dimensions[part_dimension_name]
                                new_field_definition = dict(field_definition) # Copy original field definition
                                new_field_definition['cache_variable_name'] = field_definition['short_name']
                                new_field_definition['short_name'] += '_by_' + part_dimension_name
                                new_field_definition['dimensions'] = part_dimension_name
                                new_field_definition['column_offset'] = column_offset
                                new_field_definition['columns'] = columns
                                new_field_definition['column_start_index'] += column_offset
                                new_field_definition['dimension_size'] = self.dimensions[part_dimension_name]
                                
                                field_definition_index += 1
                                logger.debug('Inserting new part-dimension field definition for {} at index {}'.format(new_field_definition['short_name'], field_definition_index))
                                self.field_definitions.insert(field_definition_index, new_field_definition)
                                column_offset += columns # Offset next column 
                            
                    else: # No default dimension name found from settings
                        match = re.match('(.+)_by_(.+)', field_name)
                        if match:
                            dimension_name = match.group(2) # This probably won't happen
                        else:
                            dimension_name = field_name + '_dim' # Nasty, but it will get the job done.  
                            
                        field_definition['dimensions'] = dimension_name
                        if not self.dimensions.get(dimension_name):
                            self.dimensions[dimension_name] = field_dimension_size
                        else:
                            assert self.dimensions[dimension_name] == field_dimension_size, 'Invalid size for variable {}. Expected {}, got {}'.format(field_name,
                                                                                                                                                       self.dimensions[dimension_name],
                                                                                                                                                       field_dimension_size
                                                                                                                                                       ) 
                        
                        logger.debug('Variable {} has dimension {} of size {}'.format(field_name,
                                                                                      dimension_name, 
                                                                                      field_dimension_size))                        
                           
                    column_start_index += field_dimension_size
                    
                else: # 1D variable
                    field_dimension_size = 0 
                    column_start_index += 1 # Single column only
                
                field_definition['dimension_size'] = field_dimension_size
                
            logger.debug('self.dimensions: {}'.format(pformat(self.dimensions)))
            logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
        
        
        def fix_all_field_precisions():
            '''
            Helper function to reduce field datatype size if there is no loss in precision
            N.B: May modify fill_value variable attribute to same format as data if fill_value is not in data
            '''
            for field_definition in self.field_definitions:
                short_name = field_definition['short_name']
                cache_variable_name = field_definition.get('cache_variable_name')
                
                # Skip dimension fields
                if short_name in self.settings['dimension_fields']:
                    continue
                
                dtype = field_definition['dtype']
                
                if cache_variable_name: # Part dimension  - shared cache variable
                    cache_variable = self._nc_cache_dataset.variables[cache_variable_name]
                    data_array = cache_variable[:,field_definition['column_offset']:field_definition['columns']]
                else: # Normal data variable
                    cache_variable = self._nc_cache_dataset.variables[short_name]
                    data_array = cache_variable[:]
                
                # Include fill value if required
                if type(data_array) == np.ma.core.MaskedArray:
                    logger.debug('Array is masked. Including fill value.')
                    fill_value = data_array.fill_value # Use masked array fill value instead of any provided value
                    no_data_mask = data_array.mask
                    data_array = data_array.data # Include mask value
                else:
                    fill_value = field_definition.get('fill_value')
                    
                if fill_value is None:
                    no_data_mask = []
                else:
                    no_data_mask = (data_array == fill_value)
            
                #logger.debug('short_name: {}, data_array: {}'.format(short_name, data_array))
                precision_change_result = fix_field_precision(data_array, dtype, 
                                                              field_definition['decimal_places'],
                                                              no_data_mask,
                                                              fill_value
                                                              ) 
                # aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format, modified_fill_value
                if precision_change_result:    
                    logger.info('Datatype for variable {} changed from {} to {}'.format(short_name, dtype, precision_change_result[1]))
                    #logger.debug('precision_change_result: {}'.format(precision_change_result))
                    field_definition['format'] = precision_change_result[0]
                    field_definition['dtype'] = precision_change_result[1]
                    field_definition['width_specifier'] = precision_change_result[3]
                    field_definition['decimal_places'] = precision_change_result[4]
                    
                    # Try truncating modified fill value
                    modified_fill_value = precision_change_result[6]
                    
                    # Update the format string in variable attributes
                    variable_attributes = field_definition.get('variable_attributes') or {}
                    if not variable_attributes:
                        field_definition['variable_attributes'] = variable_attributes
                    variable_attributes['aseg_gdf_format'] = precision_change_result[0]
                else: # Data type unchanged
                    logger.debug('Datatype for variable {} not changed'.format(short_name))
                    # Try truncating original fill value
                    modified_fill_value = truncate(fill_value, 
                                                   data_array, 
                                                   no_data_mask, 
                                                   field_definition['width_specifier'], 
                                                   field_definition['decimal_places'])
                    
                #logger.debug('(fill_value: {} modified_fill_value: {}'.format(fill_value, modified_fill_value))
                # Fill value has been modified                     
                if (fill_value is not None 
                    and modified_fill_value is not None 
                    and (fill_value != modified_fill_value)): 
                    
                    field_definition['fill_value'] = modified_fill_value
                    # Change fill_value in array variable
                    data_array[np.logical_or((data_array == fill_value), 
                                             np.isnan(data_array))] = modified_fill_value
                    cache_variable[:] = data_array
                
                    logger.info('{} fill_value changed from {} to {} for format {}'.format(short_name,
                                                                             fill_value, 
                                                                             modified_fill_value,
                                                                             field_definition['format'])
                                                                             )
            
                
                
        # Start of actual __init__() definition
        self.nc_cache_path = None
        self._nc_cache_dataset = None
        self.column_count = None # Number of columns in .dat file
        self.space_delimited = space_delimited

        if crs_string:
            self.spatial_ref = get_spatial_ref_from_wkt(crs_string)
            logger.debug('CRS set from supplied crs_string {}'.format(crs_string))
        else:
            self.spatial_ref = None

        self.dimensions = OrderedDict()

        self.settings_path = settings_path or os.path.join(os.path.dirname(__file__), 
                                                           'aseg_gdf_settings.yml')        
        try:
            self.settings = yaml.safe_load(open(self.settings_path))
        except:
            self.settings = {}
            
        logger.debug('self.settings: {}'.format(pformat(self.settings)))

        ToNetCDFConverter.__init__(self, 
                                 nc_out_path, 
                                 netcdf_format, 
                                 default_chunk_size=default_chunk_size, 
                                 default_variable_parameters=default_variable_parameters
                                 )
        
        self.aem_dat_path = aem_dat_path
        self.dfn_path = dfn_path
        
        try:
            # Parse .dfn file and apply overrides
            get_field_definitions() 
    
            # Read data from self.dfn_path into temporary netCDF variable
            read_data_file()
            
            # Needs to be done after data file has been read in order to access dimension data in first row
            modify_field_definitions()
            
            # Fix excessive precision if required - N.B: Will change field datatypes if no loss in precision
            if fix_precision:
                fix_all_field_precisions()
        except Exception as e:
            logger.error('Unable to create ASEGGDF2NetCDFConverter object: {}'.format(e))
            self.__del__()
                       
    def __del__(self):
        '''
        Destructor for class ASEGGDF2NetCDFConverter
        Closes and removes temporary cache file
        '''
        try:
            self._nc_cache_dataset.close()
            logger.debug('Closed temporary cache file {}'.format(self.nc_cache_path))
        except:
            pass
         
        try:
            os.remove(self.nc_cache_path) 
            logger.debug('Deleted temporary cache file {}'.format(self.nc_cache_path))
        except:
            pass
        
        ToNetCDFConverter.__del__(self)

    
    def get_raw_data(self, variable_name):
        '''
        Helper function to return array corresponding to short_name from self._nc_cache_dataset
        '''
        try:
            return self._nc_cache_dataset.variables[variable_name][:]
        except:
            try: # Try part dimension with shared cache variable
                field_definition = [field_definition 
                                    for field_definition in self.field_definitions
                                    if field_definition.get('short_name') == variable_name
                                    ][0]
                return self._nc_cache_dataset.variables[field_definition['cache_variable_name']][:,field_definition['column_offset']:field_definition['column_offset']+field_definition['columns']]
            except:
                return None        
        
        
    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs       
        '''
        #TODO: implement search lists for different variable names
        
        metadata_dict = {'title': 'Dataset read from ASEG-GDF file {}'.format(os.path.basename(self.aem_dat_path)),
            'Conventions': "CF-1.6,ACDD-1.3",
            'featureType': "trajectory",
            #TODO: Sort out standard names for elevation and get rid of the DTM case. Should this be min(elevation-DOI)?
            'geospatial_vertical_min': np.min(self.get_raw_data('elevation')) or np.min(self.get_raw_data('dtm')),
            'geospatial_vertical_max': np.max(self.get_raw_data('elevation')) or np.min(self.get_raw_data('dtm')), 
            'geospatial_vertical_units': "m",
            'geospatial_vertical_resolution': "point",
            'geospatial_vertical_positive': "up",
            'history': 'Converted from ASEG-GDF file {} using definitions file {}'.format(self.aem_dat_path,
                                                                                     self.dfn_path),
            'date_created': datetime.now().isoformat(),
            'geospatial_east_resolution': "point",
            'geospatial_north_resolution': "point",
            }
            
        try:
            point_count = self.nc_output_dataset.dimensions['point'].size
            coordinates = np.zeros(shape=(point_count, 2), dtype=np.float64)
            
            if set(['longitude', 'latitude']) <= set(self.nc_output_dataset.variables.keys()):
                coordinates[:,0] = self.nc_output_dataset.variables['longitude'][:]
                coordinates[:,1] = self.nc_output_dataset.variables['latitude'][:]
                
                metadata_dict.update({
                    'geospatial_lon_min': np.min(coordinates[:,0]),
                    'geospatial_lon_max': np.max(coordinates[:,0]),
                    'geospatial_lon_units': "degrees East",
                    'geospatial_lat_min': np.min(coordinates[:,1]),
                    'geospatial_lat_max': np.max(coordinates[:,1]),
                    'geospatial_lat_units': "degrees North",
                    })


            elif set(['easting', 'northing']) <= set(self.nc_output_dataset.variables.keys()): # CRS is in UTM
                coordinates[:,0] = self.nc_output_dataset.variables['easting'][:]
                coordinates[:,1] = self.nc_output_dataset.variables['northing'][:]
                
                metadata_dict.update({
                    'geospatial_east_min': np.min(coordinates[:,0]),
                    'geospatial_east_max': np.max(coordinates[:,0]),
                    'geospatial_east_units': "m",
                    'geospatial_north_min': np.min(coordinates[:,1]),
                    'geospatial_north_max': np.max(coordinates[:,1]),
                    'geospatial_north_units': "m",
                    })
            
            else:
                raise BaseException('Unrecognised coordinates')
            
            #Compute convex hull and add GML representation to metadata
            #logger.debug('coordinates: {}'.format(coordinates))
            convex_hull = points2convex_hull(coordinates)        
            metadata_dict['geospatial_bounds'] = 'POLYGON((' + ', '.join([' '.join(
                ['%.4f' % ordinate for ordinate in coordinates]) for coordinates in convex_hull]) + '))'

        except BaseException as e:
            logger.warning('Unable to set spatial attributes: {}'.format(e))
            
        if self.settings.get('keywords'):
            metadata_dict['keywords'] = self.settings['keywords']

        logger.debug('metadata_dict: {}'.format(metadata_dict))
        return metadata_dict
    
    
    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs       
        '''
        return self.dimensions
    
    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects       
        ''' 
        def index_variable_generator():
            '''
            Helper generator to yield indexing variables
            N.B: Has side effect of creating one new dimension for each index field. 
            This new dimension will NOT appear in self.dimensions
            '''
            # Process index variables
            try:
                logger.info('{} {} values found'.format(len(lookup_array), field_definition['short_name']))
                #logger.debug('lookup_array: {},\nindex_start_indices: {},\nindex_point_counts: {}'.format(lookup_array, index_start_indices, index_point_counts))
            except:
                logger.info('Unable to create {} indexing variables'.format(field_definition['short_name']))
                return   

            # This is a slightly ugly side effect
            logger.info('\tCreating dimension for {}'.format(field_definition['short_name']))
            self.nc_output_dataset.createDimension(dimname=field_definition['short_name'], 
                                                   size=len(lookup_array))
            
            logger.info('\tWriting {} indexing variables'.format(field_definition['short_name']))
            
            logger.info('\t\tWriting {} values'.format(field_definition['short_name']))
            
            variable_attributes = (field_definition.get('variable_attributes') or {})
            if field_definition.get('long_name'):
                variable_attributes['long_name'] = field_definition['long_name'] 
                
            yield NetCDFVariable(short_name=field_definition['short_name'], 
                                 data=lookup_array, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes=variable_attributes, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting index of first point in each {}'.format(field_definition['short_name']))
            yield NetCDFVariable(short_name='{}_start_index'.format(field_definition['short_name']), 
                                 data=index_start_indices, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'zero-based index of the first point in each {}'.format(field_definition['short_name'])}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting point count for each {}'.format(field_definition['short_name']))
            yield NetCDFVariable(short_name='{}_point_count'.format(field_definition['short_name']), 
                                 data=index_point_counts, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'number of points in each {}'.format(field_definition['short_name'])}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
              
        def lookup_variable_generator():
            '''
            Helper generator to yield lookup variables for string fields
            N.B: Has side effect of creating one new dimension for each lookup field. 
            This new dimension will NOT appear in self.dimensions
            '''
            assert lookup_array.shape[0] > 1, 'lookup_array must have more than one value'
            
            short_name = field_definition['short_name']

            logger.info('\tWriting {} lookup variables'.format(short_name))
            
            variable_attributes = (field_definition.get('variable_attributes') or {})
            if field_definition.get('long_name'):
                variable_attributes['long_name'] = field_definition['long_name'] 
                
            # Creating a new dimension is a slightly ugly side effect
            logger.info('\tCreating dimension for {}'.format(short_name))
            self.nc_output_dataset.createDimension(dimname=short_name, 
                                                   size=len(lookup_array))
            
            logger.info('\t\tWriting {} {} lookup values to array variable {}'.format(len(lookup_array), 
                                                                                      short_name, 
                                                                                      short_name))
            yield NetCDFVariable(short_name=short_name, 
                                 data=lookup_array, 
                                 dimensions=[short_name], 
                                 fill_value=None, 
                                 attributes=variable_attributes, 
                                 dtype='str' if (lookup_array.dtype == object) else lookup_array.dtype,
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            variable_name = '{}_index'.format(short_name)
            logger.info('\t\tWriting {} lookup indices to array variable {}'.format(short_name,
                                                                                    variable_name))
            yield NetCDFVariable(short_name=variable_name, 
                                 data=index_array, 
                                 dimensions=['point'], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'zero-based index of value in {}'.format(short_name),
                                             'lookup': short_name
                                             }, 
                                 dtype='int8' if len(lookup_array) < 128 else 'int32' if len(lookup_array) < 32768 else index_array.dtype,
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )

            return        
        
        
        def get_scalar_variable():
            '''
            Helper function to return scalar variable for single-value field
            '''
            assert lookup_array.shape[0] == 1, 'lookup_array must have exactly one (i.e. unique) value'
            short_name = field_definition['short_name']

            logger.info('\tWriting single {} value to scalar variable'.format(short_name))
            
            variable_attributes = (field_definition.get('variable_attributes') or {})
            if field_definition.get('long_name'):
                variable_attributes['long_name'] = field_definition['long_name'] 
                
            return NetCDFVariable(short_name=short_name, 
                                 data=lookup_array, 
                                 dimensions=[], 
                                 fill_value=None, 
                                 attributes=variable_attributes, 
                                 dtype='str' if lookup_array.dtype == object else lookup_array.dtype,
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
            
                        
        # Start of variable_generator function
        # Create crs variable
        yield self.build_crs_variable(self.spatial_ref)
        
        # Create and yield variables
        for field_definition in self.field_definitions:
            #logger.debug('field_definition: {}'.format(pformat(field_definition)))
            
            short_name = field_definition['short_name']
            
            # Skip dimension fields or ignored fields
            if field_definition['short_name'] in self.settings['dimension_fields'] + self.settings['ignored_fields']:
                logger.debug('\tIgnoring variable {}'.format(short_name))
                continue
            
            field_attributes = OrderedDict()
            
            raw_data = self.get_raw_data(short_name)
            
            lookup_array, index_start_indices, index_array, index_point_counts = np.unique(raw_data, 
                                                                              return_index=True, 
                                                                              return_inverse=True, 
                                                                              return_counts=True
                                                                              )            
            
            unique_value_count = lookup_array.shape[0]
            logger.debug('{} unique values found for {}'.format(unique_value_count, short_name))

            # Only one unique value in data for ALL points in a 1D variable - write to scalar variable
            #TODO: Check what happens with masked arrays
            if (field_definition['columns'] == 1) and (unique_value_count == 1) and (index_point_counts[0] == raw_data.shape[0]): 
                yield get_scalar_variable()
                continue
            
            # Process string fields or designated lookup fields as lookups
            elif (field_definition['format'].startswith('A') # String field
                or short_name in self.settings['lookup_fields'] # Designated lookup field
                ):
                for lookup_variable in lookup_variable_generator():
                    # Create index dimension
                    yield lookup_variable
                     
                continue
            
            #===================================================================
            # # THIS HAS BEEN REMOVED BECAUSE WE NO LONGER USE INDEXING VARIABLES
            # # Process designated index fields
            # elif short_name in self.settings['index_fields']:
            #     for index_variable in index_variable_generator():
            #         # Create index dimension
            #         yield index_variable
            #         
            #     continue
            #===================================================================
            
            # Process normal data field
            long_name = field_definition.get('long_name')
            if long_name:
                field_attributes['long_name'] = long_name
                
            units = field_definition.get('units')
            if units:
                field_attributes['units'] = units
                
            # Append any recognised variable attributes read from .dfn
            field_attributes.update(field_definition.get('variable_attributes') or {})
                
            dtype = field_definition.get('dtype')
                
            fill_value = field_definition.get('fill_value')
                
            if not field_definition['dimension_size']: # 1D Variable
                logger.info('\tWriting 1D {} variable {}'.format(dtype, short_name))
                
                yield NetCDFVariable(short_name=short_name, 
                                     data=self.get_raw_data(short_name), 
                                     dimensions=['point'], 
                                     fill_value=fill_value, 
                                     attributes=field_attributes, 
                                     dtype=dtype,
                                     chunk_size=self.default_chunk_size,
                                     variable_parameters=self.default_variable_parameters
                                     )
        
            #=======================================================================
            # # Create bad_data_mask array from depth of investigation
            # top_depth = self.get_raw_data('layer_top_depth')
            # depth_of_investigation = self.get_raw_data('depth_of_investigation')
            # bad_data_mask = top_depth > np.repeat(depth_of_investigation[:, np.newaxis], 
            #                                     top_depth.shape[1], 
            #                                     axis=1)
            # logger.debug('{} bad conductivity values found for masking'.format(np.count_nonzero(bad_data_mask)))
            #=======================================================================
            else: # 2D variable
                #TODO: Move this code to post-processing
                # Convert resistivity to conductivity
                if short_name == 'resistivity':
                    short_name = 'conductivity'
                    field_attributes = {'long_name': 'Layer conductivity', 'units': 'S/m'}
                    data_array =  1.0 / self.get_raw_data('resistivity')
                    fill_value = 0
                    #===========================================================
                    # data_array[bad_data_mask] = fill_value
                    #===========================================================
                    logger.debug('\nconductivity: {}'.format(data_array))
                elif short_name == 'resistivity_uncertainty':
                    # Convert resistivity_uncertainty to absolute_conductivity_uncertainty
                    #TODO: Check whether resistivity_uncertainty is a proportion or a percentage - the former is assumed
                    # Search for "reciprocal" in http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf
                    short_name = 'conductivity_uncertainty'
                    field_attributes = {'long_name': 'Absolute uncertainty of layer conductivity', 'units': 'S/m'}
                    data_array =  self.get_raw_data('resistivity_uncertainty') / self.get_raw_data('resistivity')
                    fill_value = 0
                    #===========================================================
                    # data_array[bad_data_mask] = fill_value
                    #===========================================================
    
                    #===========================================================
                    # logger.debug('\nresistivity_uncertainty: {}'.format(self.get_raw_data('resistivity_uncertainty')))
                    # logger.debug('\nresistivity: {}'.format(self.get_raw_data('resistivity')))
                    # logger.debug('\nabsolute_resistivity_uncertainty: {}'.format(self.get_raw_data('resistivity') 
                    #                                                            * self.get_raw_data('resistivity_uncertainty')))
                    # logger.debug('\nabsolute_conductivity_uncertainty: {}\n'.format(data_array))
                    #===========================================================
                else:
                    # Don't mess with values
                    data_array = self.get_raw_data(short_name)
                    fill_value = None
                                          
                logger.info('\tWriting 2D {} variable {}'.format(dtype, short_name))
    
                yield NetCDFVariable(short_name=short_name, 
                                     data=data_array, 
                                     dimensions=['point', field_definition['dimensions']], 
                                     fill_value=fill_value, 
                                     attributes=field_attributes, 
                                     dtype=dtype,
                                     chunk_size=self.default_chunk_size,
                                     variable_parameters=self.default_variable_parameters
                                     )
            
        return
    
    def postprocess_netcdf(self):
        '''
        Function to perform post-processing on netCDF file after dimensions and variables
        have been created. Overrides base class.
        
        This will create crs, longitude and latitude variables for unprojected CRS, 
        and recompute global attributes
        '''
        default_crs_wkt = self.settings['default_crs_wkt']
        
        crs_var = self.nc_output_dataset.variables.get('crs')
        transverse_mercator_var = self.nc_output_dataset.variables.get('transverse_mercator')
        logger.debug('crs_var: {}'.format(crs_var))
        logger.debug('transverse_mercator_var: {}'.format(transverse_mercator_var))
        
        # If dataset has UTM coordinates and not unprojected ones
        if (transverse_mercator_var is not None 
            and crs_var is None
            and set(['easting', 'northing']) <= set(self.nc_output_dataset.variables.keys())
            ):
            # Build GDA94 crs variable and write it to self.nc_output_dataset
            logger.info('Creating crs, longitude and latitude variables for unprojected CRS')
            point_count = self.nc_output_dataset.dimensions['point'].size
            logger.debug('point_count: {}'.format(point_count))
            
            utm_coords = np.ones(shape=(point_count, 2), dtype=np.float64) * -999
            utm_coords[:,0] = self.nc_output_dataset.variables['easting'][:]
            utm_coords[:,1] = self.nc_output_dataset.variables['northing'][:]
            
            gda94_coords = transform_coords(utm_coords, 
                                            transverse_mercator_var.spatial_ref, 
                                            default_crs_wkt
                                            )
            
            # Create and write crs variable
            logger.info('Creating new crs variable for unprojected CRS')
            self.build_crs_variable(get_spatial_ref_from_wkt(default_crs_wkt)
                                    ).create_var_in_dataset(self.nc_output_dataset)
                                    
            # Create and write longitude variable
            logger.info('Creating new longitude variable')
            NetCDFVariable('longitude', 
                           gda94_coords[:,0], 
                           ['point'], 
                           fill_value=-999, 
                           attributes={'long_name': 'Longitude', 'units': 'degrees East'}
                           ).create_var_in_dataset(self.nc_output_dataset)
            
            # Create and write latitude variable
            logger.info('Creating new latitude variable')
            NetCDFVariable('latitude', 
                           gda94_coords[:,1], 
                           ['point'], 
                           fill_value=-999, 
                           attributes={'long_name': 'Latitude', 'units': 'degrees North'}
                           ).create_var_in_dataset(self.nc_output_dataset)
                           
            logger.info('Re-writing new global attributes for new CRS')
            for attribute_name, attribute_value in iter(self.get_global_attributes().items()):
                setattr(self.nc_output_dataset, attribute_name, attribute_value or '')
            
    

def main():
    '''
    Main function
    '''
    def get_args():
        """
        Handles all the arguments that are passed into the script

        :return: Returns a parsed version of the arguments.
        """
        parser = argparse.ArgumentParser(description='Convert ASEG-GDF file to netCDF')
        parser.add_argument("-f", "--dfn",
                            help="Path to .dfn file",
                            type=str,
                            dest="dfn_in_path")
        parser.add_argument("-s", "--settings",
                            help="Path to settings file",
                            type=str,
                            dest="settings_path")
        parser.add_argument("-r", "--crs",
                            help="Coordinate Reference System string (e.g. GDA94, EPSG:4283)",
                            type=str,
                            dest="crs")
        parser.add_argument("-c", "--chunking",
                            help="Chunking size in each dimension",
                            type=int,
                            dest="chunking",
                            default=1024)
        parser.add_argument('-o', '--optimise', 
                            help='Optimise datatypes to reduce variable size without affecting precision. Default is True',
                            type=int,
                            dest='optimise',
                            default=True
                            )
        parser.add_argument('-a', '--space_delimited', 
                            help='Read .dat file as space-delimited instead of fixed-length. Default is False',
                            type=int,
                            dest='space_delimited',
                            default=False
                            )
        
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='output debug information. Default is no debug info')
        
        parser.add_argument('positional_args', 
                            nargs=argparse.REMAINDER,
                            help='<dat_in_path> [<nc_out_path>]')

        return parser.parse_args()
    
    args = get_args()
    
    # Setup Logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    base_logger = logging.getLogger('geophys_utils.netcdf_converter') # Logger for base class ToNetCDFConverter
    base_logger.setLevel(level=log_level)
    logger.setLevel(level=log_level)

    assert 1 <= len(args.positional_args) <= 2, 'Invalid number of positional arguments.\n\
Usage: python {} <options> <dat_in_path> [<nc_out_path>]'.format(os.path.basename(sys.argv[0]))

    dat_in_path = args.positional_args[0] # 'C:\\Temp\\Groundwater Data\\ord_bonaparte_nbc_main_aquifer_clipped.dat'

    if len(args.positional_args) >= 2:
        nc_out_path = args.positional_args[1] # 'C:\\Temp\\dat_test.nc'
    else:
        nc_out_path = os.path.splitext(dat_in_path)[0] + '.nc'
        
    dfn_in_path = args.dfn_in_path or os.path.splitext(dat_in_path)[0] + '.dfn'
       
    try:
        d2n = ASEGGDF2NetCDFConverter(nc_out_path, 
                                      dat_in_path, 
                                      dfn_in_path, 
                                      crs_string=args.crs,
                                      default_chunk_size=args.chunking,
                                      settings_path=args.settings_path,
                                      fix_precision=args.optimise,
                                      space_delimited=args.space_delimited
                                      )
        d2n.convert2netcdf()
        logger.info('Finished writing netCDF file {}'.format(nc_out_path))
    except Exception as e:
        logger.error('Failed to create netCDF file {}'.format(e))
        try:
            del d2n
        except Exception as e:
            logger.debug('Unable to delete ASEGGDF2NetCDFConverter object: {}'.format(e))
        
        try:
            os.remove(nc_out_path) 
            logger.debug('Deleted incomplete netCDF file {}'.format(nc_out_path))
        except Exception as e:
            logger.debug('Unable to delete incomplete netCDF file {}: {}'.format(nc_out_path, e))
            
        
    
    #===========================================================================
    # logger.debug('\nGlobal attributes:')
    # logger.debug(pformat(d2n.nc_output_dataset.__dict__))
    # logger.debug('Dimensions:')
    # logger.debug(pformat(d2n.nc_output_dataset.dimensions))
    # logger.debug('Variables:')
    # for variable_name in d2n.nc_output_dataset.variables.keys():
    #     variable = d2n.nc_output_dataset.variables[variable_name]
    #     logger.debug(pformat(variable))
    #     logger.debug(variable[:])
    #===========================================================================

if __name__ == '__main__':
    # Setup logging handlers if required
    if not logger.handlers:
        # Set handler for root logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        #console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for logger {}'.format(logger.name))

    main()
