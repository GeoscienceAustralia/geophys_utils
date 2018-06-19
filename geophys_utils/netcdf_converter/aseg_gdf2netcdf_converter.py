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
from geophys_utils.netcdf_converter.aseg_gdf_format_dtype import aseg_gdf_format2dtype

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

TEMP_DIR = tempfile.gettempdir()
#TEMP_DIR = 'E:\Temp'

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
                 netcdf_format='NETCDF4_CLASSIC', 
                 default_chunk_size=None, 
                 default_variable_parameters=None,
                 settings_path=None
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
        '''

        def get_field_definitions():
            '''
            Function to read raw field definitions from .dfn file
            Will set self.dimensions as an Ordereddict of dimension sizes keyed by dimension name
            '''
            def parse_dfn_file(dfn_path):
                logger.info('Reading definitions file {}'.format(dfn_path))
                field_definitions = []
                crs=None
                
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
    
                    #logger.debug('key_value_pairs: {},\npositional_value_list: {}'.format(pformat(key_value_pairs), pformat(positional_value_list))) 
                
                    # Column definition
                    if key_value_pairs.get('RT') == '' and (positional_value_list 
                                                            and positional_value_list[0] != 'END DEFN'): 
                        short_name = positional_value_list[0].lower()
                        fmt = positional_value_list[1] if len(positional_value_list) >= 2 else None
                        units = key_value_pairs.get('UNITS') or key_value_pairs.get('UNIT')
                        long_name = key_value_pairs.get('NAME') or (positional_value_list[2] if len(positional_value_list) >= 3 else None)
                        fill_value = key_value_pairs.get('NULL')
                        
                        # Parse format to determine columns, data type and numeric format
                        dtype, columns, integer_digits, fractional_digits = aseg_gdf_format2dtype(fmt)
                                    
                        field_dict = {'short_name': short_name,
                                      'format': fmt,
                                      'long_name': long_name,
                                      'dtype': dtype,
                                      'columns': columns,
                                      'integer_digits': integer_digits, 
                                      'fractional_digits': fractional_digits
                                      }
                        if units:
                            field_dict['units'] = units
                        if fill_value is not None:
                            field_dict['fill_value'] = fill_value
                             
                        field_definitions.append(field_dict)    
                                    
                    # Projection definition
                    elif key_value_pairs.get('RT') == 'PROJ': 
                        #logger.debug('split_lists[1]: {}'.format(pformat(split_lists[1])))
                         
                        #TODO: Parse projection properly
                        if (key_value_pairs.get('PROJNAME') == 'GDA94 / MGA zone 54'):
                             
                            # TODO: Remove this hard-coded hack
                            wkt = 'EPSG:28354'
                            crs = get_spatial_ref_from_wkt(wkt)
                            break # Nothing more to do
                
                    if not crs: #TODO: Fix this hard-coded hack
                        crs = get_spatial_ref_from_wkt('EPSG:28353') # GDA94 / MGA zone 53
                
                return field_definitions, crs
            
            self.field_definitions, self.crs = parse_dfn_file(dfn_path)
            
            # Read overriding field definition values from settings
            if self.settings.get('field_definitions'):
                for field_definition in self.field_definitions:
                    overriding_field_definition = self.settings['field_definitions'].get(field_definition['short_name'])
                    if overriding_field_definition:
                        field_definition.update(overriding_field_definition)
            
            logger.debug('self.dimensions: {}'.format(pformat(self.dimensions)))
            logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))

        
        def read_data_file(): 
            '''
            Function to read data file into temporary netCDF cache
            Converts entire numeric file into a single variable of float32 datatype
            This is done to permit efficient column-wise data access but it will fail for any
            invalid floating point values in data file (except in possible header row which is ignored)
            '''   
            def create_nc_cache():
                '''
                Function to create temporary cache file with one 2D variable
                Needs to have self.column_count defined to work
                '''
                assert self.column_count, 'self.column_count not defined'
                
                self.nc_cache_path = os.path.join(TEMP_DIR, re.sub('\W+', '_', os.path.splitext(self.aem_dat_path)[0]) + '.nc')
                self._nc_cache_dataset = netCDF4.Dataset(self.nc_cache_path, mode="w", clobber=True, format='NETCDF4')
                self._nc_cache_dataset.createDimension(dimname='columns', size=self.column_count)
                self._nc_cache_dataset.createDimension(dimname='rows', size=None)
                self.cache_variable = self._nc_cache_dataset.createVariable('cache',
                                                                            'float64',
                                                                            ('rows', 'columns'),
                                                                            chunksizes=(CACHE_CHUNK_ROWS, 1), # Chunk by column
                                                                            **NetCDFVariable.DEFAULT_VARIABLE_PARAMETERS
                                                                            )
                
                logger.debug('Created temporary cache file {}'.format(self.nc_cache_path))
            
            logger.info('Reading data file {}'.format(aem_dat_path))
            aem_dat_file = open(aem_dat_path, 'r')
             
            self.column_count = 0
            self.total_points = 0
            
            for line in aem_dat_file:
                line = line.strip()
                if not line: # Skip empty lines
                    continue
                
                try:
                    row_list = [float(value) for value in re.sub('\s+', ',', line).split(',')]
                except ValueError:
                    # Ignore potential non-numeric header row, otherwise raise exception for non-numeric data
                    if self.total_points:
                        raise
                    else:
                        continue
                #logger.debug('row: {}'.format(row))
                
                if self.column_count:
                    assert self.column_count == len(row_list), 'Inconsistent column count. Expected {}, found {}.'.format(self.column_count, 
                                                                                                                   len(row_list)
                                                                                                                   )
                else: # First row
                    self.column_count = len(row_list)
                    logger.debug('self.column_count: {}'.format(self.column_count))
                    create_nc_cache() # Create temporary netCDF file with appropriately dimensioned variable
                    memory_cache_array = np.zeros(shape=(CACHE_CHUNK_ROWS, self.column_count), dtype='float64')
                
                try:
                    # Write row to memory cache
                    memory_cache_array[self.total_points%CACHE_CHUNK_ROWS,:] = np.array(row_list, dtype='float64')
                except Exception as e:
                    logger.debug('Error with row = {}'.format(row_list))
                    raise e
                
                self.total_points += 1
                
                if self.total_points % CACHE_CHUNK_ROWS == 0:
                    # Write memory cache to disk cache
                    logger.debug('Writing rows {}-{} to disk cache'.format(self.total_points-CACHE_CHUNK_ROWS, self.total_points-1))
                    self.cache_variable[self.total_points-CACHE_CHUNK_ROWS:self.total_points,:] = memory_cache_array
                
                if self.total_points == self.total_points // 10000 * 10000:
                    logger.info('{} points read'.format(self.total_points))
                    #logger.debug('line: {}'.format(line))
                    
                if POINT_LIMIT and self.total_points >= POINT_LIMIT:
                    logger.debug('Truncating input for testing after {} points'.format(POINT_LIMIT))
                    break
                
            logger.debug('Writing final rows {}-{} to disk cache'.format(self.total_points-self.total_points%CACHE_CHUNK_ROWS, self.total_points-1))
            self.cache_variable[self.total_points-self.total_points%CACHE_CHUNK_ROWS:self.total_points,:] = memory_cache_array[0:self.total_points%CACHE_CHUNK_ROWS,:]
            
            aem_dat_file.close()             
    
            logger.info('A total of {} points were read'.format(self.total_points))
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
                logger.debug('field_definition_index: {}, field_name: {}'.format(field_definition_index, field_name))
                
                if field_name in self.settings['dimension_fields']:
                    self.dimensions[field_name] = int(self.cache_variable[0, column_start_index]) # Read value from first line
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
                    logger.debug('default_dimension_names: {}'.format(default_dimension_names))
                    if default_dimension_names: # Default dimension name(s) found from settings
                        if type(default_dimension_names) == str: # Only one dimension specified as a string
                            actual_dimension_size = self.dimensions.get(default_dimension_names)
                            logger.debug('actual_dimension_size: {}'.format(actual_dimension_size))
                            if actual_dimension_size is None: # Default dimension not already defined - create it
                                self.dimensions[default_dimension_names] = field_dimension_size
                            else:
                                assert field_dimension_size == actual_dimension_size, 'Mismatched dimension sizes'
                                
                            logger.debug('Variable {} has dimension {} of size {}'.format(field_name,
                                                                                      default_dimension_names, 
                                                                                      field_dimension_size))
                        
                        # Use dimension size values read from data file to create multiple part-arrays from a single composite
                        # N.B: This could be considered to be encouraging bad behaviour - the .dfn file should 
                        # ideally be re-written
                        if type(default_dimension_names) == list: # Multiple part-dimensions specified
                            assert set(default_dimension_names) <= set(self.dimensions.keys()), 'Dimensions {} not all defined in {}'.format(default_dimension_names, self.dimensions.keys())
                            assert sum([self.dimensions[key] for key in default_dimension_names]) == field_dimension_size, 'Incorrect dimension sizes'
                            logger.debug('Removing original composite field definition')
                            self.field_definitions.remove(field_definition)
                            field_definition_index -= 1
                            # Create new field definitions for part dimensions
                            column_offset = 0
                            for part_dimension_name in default_dimension_names:
                                new_field_definition = dict(field_definition) # Copy original field definition
                                new_field_definition['short_name'] += '_x_' + part_dimension_name
                                new_field_definition['dimensions'] = part_dimension_name
                                new_field_definition['column_start_index'] += column_offset
                                new_field_definition['dimension_size'] = self.dimensions[part_dimension_name]
                                
                                field_definition_index += 1
                                logger.debug('Inserting new part-dimension field definition for {} at index {}'.format(new_field_definition['short_name'], field_definition_index))
                                self.field_definitions.insert(field_definition_index, new_field_definition)
                                column_offset += self.dimensions[part_dimension_name] # Offset next column 
                            
                    else: # No default dimension name found from settings
                        match = re.match('(.+)_x_(.+)', field_name)
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
        
            
        # Start of actual __init__() definition
        self.nc_cache_path = None
        self._nc_cache_dataset = None
        self.cache_variable = None
        self.column_count = None # Number of columns in .dat file

        self.dimensions = OrderedDict()

        self.settings_path = settings_path or os.path.splitext(__file__)[0] + '_settings.yml'
        
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
        
        # Parse .dfn file and apply overrides
        get_field_definitions() 

        # Read data from self.dfn_path into temporary netCDF variable
        read_data_file()
        
        # Needs to be done after data file has been read in order to access dimension data in first row
        modify_field_definitions()
        
                       
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

    
    def get_raw_data(self, short_name):
        '''
        Helper function to return array corresponding to short_name from self.cache_variable
        '''
        try:
            field_definition = [field_definition 
                                for field_definition in self.field_definitions 
                                if field_definition['short_name'] == short_name
                                ][0]
             
            logger.debug('field_definition: {}'.format(pformat(field_definition)))
            column_start_index = field_definition['column_start_index']
            column_end_index = column_start_index + max(1, field_definition['dimension_size'])
            
            return self.cache_variable[:,column_start_index:column_end_index]
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
            'geospatial_east_min': np.min(self.get_raw_data('easting')) or np.min(self.get_raw_data('longitude')),
            'geospatial_east_max': np.max(self.get_raw_data('easting')) or np.max(self.get_raw_data('longitude')),
            'geospatial_east_units': "m",
            'geospatial_east_resolution': "point",
            'geospatial_north_min': np.min(self.get_raw_data('northing')) or np.min(self.get_raw_data('latitude')),
            'geospatial_north_max': np.max(self.get_raw_data('northing')) or np.min(self.get_raw_data('latitude')),
            'geospatial_north_units': "m",
            'geospatial_north_resolution': "point",
            #TODO: Sort out standard names for elevation and get rid of the DTM case. Should this be min(elevation-DOI)?
            'geospatial_vertical_min': np.min(self.get_raw_data('elevation')) or np.min(self.get_raw_data('DTM')),
            'geospatial_vertical_max': np.max(self.get_raw_data('elevation')) or np.min(self.get_raw_data('DTM')), 
            'geospatial_vertical_units': "m",
            'geospatial_vertical_resolution': "point",
            'geospatial_vertical_positive': "up",
            'history': 'Converted from ASEG-GDF file {} using definitions file {}'.format(self.aem_dat_path,
                                                                                     self.dfn_path),
            'date_created': datetime.now().isoformat()
            }
        
        if self.settings.get('keywords'):
            metadata_dict['keywords'] = self.settings['keywords']

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
        def index_variable_generator(field_definition):
            '''
            Helper generator to yield indexing variables
            N.B: Has side effect of creating one new dimension for each index field. 
            This new dimension will NOT appear in self.dimensions
            '''
            # Process index variables
            try:
                index_data = self.get_raw_data(field_definition['short_name'])
                
                index_values, index_start_indices, index_point_counts = np.unique(index_data, 
                                                                                  return_index=True, 
                                                                                  return_inverse=False, 
                                                                                  return_counts=True
                                                                                  )
                logger.info('{} {} values found'.format(len(index_values), field_definition['short_name']))
                logger.debug('index_values: {},\nindex_start_indices: {},\nindex_point_counts: {}'.format(index_values, index_start_indices, index_point_counts))
            except:
                logger.info('Unable to create {} indexing variables'.format(field_definition['short_name']))
                return   

            # This is a slightly ugly side effect
            logger.info('\tCreating dimension for {}'.format(field_definition['short_name']))
            self.nc_output_dataset.createDimension(dimname=field_definition['short_name'], 
                                                   size=len(index_values))
            
            logger.info('\tWriting {} indexing variables'.format(field_definition['short_name']))
            
            logger.info('\t\tWriting {} values'.format(field_definition['short_name']))
            yield NetCDFVariable(short_name=field_definition['short_name'], 
                                 data=index_values, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes={'long_name': field_definition['long_name']} if field_definition.get('long_name') else {}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting index of first point in each {}'.format(field_definition['short_name']))
            yield NetCDFVariable(short_name='index_' + field_definition['short_name'], 
                                 data=index_start_indices, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'zero-based index of the first point in each {}'.format(field_definition['short_name'])}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting point count for each {}'.format(field_definition['short_name']))
            yield NetCDFVariable(short_name='index_count_' + field_definition['short_name'], 
                                 data=index_point_counts, 
                                 dimensions=[field_definition['short_name']], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'number of points in each {}'.format(field_definition['short_name'])}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
              
        # Create crs variable
        yield self.build_crs_variable(self.crs)
        
        # Create and yield variables
        for field_definition in self.field_definitions:
            short_name = field_definition['short_name']
            
            # Skip dimension fields or ignored fields
            if field_definition['short_name'] in self.settings['dimension_fields'] + self.settings['ignored_fields']:
                logger.debug('\tIgnoring variable {}'.format(short_name))
                continue
            
            field_attributes = {}
            
            # Process index fields
            if short_name in self.settings['index_fields']:
                for index_variable in index_variable_generator(field_definition):
                    # Create index dimension
                    yield index_variable
                    
                continue
            
            long_name = field_definition.get('long_name')
            if long_name:
                field_attributes['long_name'] = long_name
                
            units = field_definition.get('units')
            if units:
                field_attributes['units'] = units
                
            dtype = field_definition.get('dtype')
                
            fill_value = field_definition.get('fill_value')
                
            logger.debug('field_definition: {}'.format(pformat(field_definition)))
            
            if not field_definition['dimension_size']: # 1D Variable
                logger.info('\tWriting 1D {} variable {}'.format(dtype, short_name))
                
                yield NetCDFVariable(short_name=short_name, 
                                     data=self.cache_variable[:,field_definition['column_start_index']], 
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
    
                    logger.debug('\nresistivity_uncertainty: {}'.format(self.get_raw_data('resistivity_uncertainty')))
                    logger.debug('\nresistivity: {}'.format(self.get_raw_data('resistivity')))
                    logger.debug('\nabsolute_resistivity_uncertainty: {}'.format(self.get_raw_data('resistivity') 
                                                                               * self.get_raw_data('resistivity_uncertainty')))
                    logger.debug('\nabsolute_conductivity_uncertainty: {}\n'.format(data_array))
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
        '''
        return
    

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
        parser.add_argument("-c", "--chunking",
                            help="Chunking size in each dimension",
                            type=int,
                            dest="chunking",
                            default=1024)
        
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
       
    d2n = ASEGGDF2NetCDFConverter(nc_out_path, 
                                  dat_in_path, 
                                  dfn_in_path, 
                                  default_chunk_size=args.chunking,
                                  settings_path=args.settings_path)
    d2n.convert2netcdf()
    logger.info('Finished writing netCDF file {}'.format(nc_out_path))
    
    logger.debug('\nGlobal attributes:')
    logger.debug(pformat(d2n.nc_output_dataset.__dict__))
    logger.debug('Dimensions:')
    logger.debug(pformat(d2n.nc_output_dataset.dimensions))
    logger.debug('Variables:')
    for variable_name in d2n.nc_output_dataset.variables.keys():
        variable = d2n.nc_output_dataset.variables[variable_name]
        logger.debug(pformat(variable))
        logger.debug(variable[:])

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
