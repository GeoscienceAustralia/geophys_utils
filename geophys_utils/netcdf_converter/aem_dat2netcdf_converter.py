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
AEMDAT2NetCDFConverter concrete class for converting AEM .dat data to netCDF

Created on 28Mar.2018

@author: Alex Ip
'''
from collections import OrderedDict
import numpy as np
import re
import os
import sys
from datetime import datetime
from pprint import pformat
import yaml
import logging

from geophys_utils.netcdf_converter import NetCDFConverter, NetCDFVariable
from geophys_utils import get_spatial_ref_from_wkt

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # Logging level for this module


class AEMDAT2NetCDFConverter(NetCDFConverter):
    '''
    AEMDAT2NetCDFConverter concrete class for converting AEM DAT data to netCDF
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
        Concrete constructor for subclass AEMDAT2NetCDFConverter
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
        #TODO: Make this a bit easier to work with - it's a bit opaque at the moment
        def parse_dfn_file(dfn_path):
            logger.info('Reading definitions file {}'.format(dfn_path))
            field_definitions = []
            crs=None
            
            dfn_file = open(dfn_path, 'r')
            for line in dfn_file:
                # generate nested lists split by ";", ":", "," then "="
                split_lists = [[[[equals_split_string.strip() for equals_split_string in comma_split_string.split('=')]
                    for comma_split_string in [comma_split_string.strip() for comma_split_string in colon_split_string.split(',')]
                    ]
                        for colon_split_string in [colon_split_string.strip() for colon_split_string in semicolon_split_string.split(':')]
                        ]
                            for semicolon_split_string in [semicolon_split_string.strip() for semicolon_split_string in line.split(';')]
                            ] 
                #logger.debug('split_lists: {}'.format(pformat(split_lists))) 
                
                try:
                    rt = split_lists[0][0][1][1] # RecordType?
                     
                    if rt == '' and len(split_lists[1]) > 1:
                        units = None
                        long_name= None
                        
                        short_name = split_lists[1][0][0][0]
                        
                        fmt = split_lists[1][1][0][0]
                        
                        if len(split_lists[1]) > 2:
                            if len(split_lists[1][2]) == 1: # Only short_name, format and long_name defined
                                long_name = split_lists[1][2][0][0]
                                
                            
                            elif len(split_lists[1][2][0]) == 2 and split_lists[1][2][0][0] == 'UNITS':
                                units = split_lists[1][2][0][1]
                                long_name = split_lists[1][2][1][0]
                            else:
                                long_name = split_lists[1][2][0][0]
                                
                        field_dict = {'short_name': short_name,
                                      'format': fmt,
                                      'long_name': long_name,
                                      }
                        if units:
                            field_dict['units'] = units
                            
                        field_definitions.append(field_dict)
                        
                    elif rt == 'PROJ':
                        #logger.debug('split_lists[1]: {}'.format(pformat(split_lists[1])))
                        
                        #TODO: Parse projection
                        if (split_lists[1][0][0][0] == 'PROJNAME' 
                            and split_lists[1][2][0][1] == 'GDA94 / MGA zone 54'):
                            
                            # TODO: Remove this hard-coded hack
                            wkt = '''PROJCS["GDA94 / MGA zone 54",
    GEOGCS["GDA94",
        DATUM["Geocentric_Datum_of_Australia_1994",
            SPHEROID["GRS 1980",6378137,298.257222101,
                AUTHORITY["EPSG","7019"]],
            TOWGS84[0,0,0,0,0,0,0],
            AUTHORITY["EPSG","6283"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4283"]],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    PROJECTION["Transverse_Mercator"],
    PARAMETER["latitude_of_origin",0],
    PARAMETER["central_meridian",141],
    PARAMETER["scale_factor",0.9996],
    PARAMETER["false_easting",500000],
    PARAMETER["false_northing",10000000],
    AUTHORITY["EPSG","28354"],
    AXIS["Easting",EAST],
    AXIS["Northing",NORTH]]'''
                            crs = get_spatial_ref_from_wkt(wkt)
                            break # Nothing more to do
                    
                except IndexError:
                    continue
            
            dfn_file.close()
            
            
            return field_definitions, crs
            
        # Start of __init__() definition
        self.settings_path = settings_path or os.path.splitext(__file__)[0] + '_settings.yml'
        
        try:
            self.settings = yaml.safe_load(open(self.settings_path))
        except:
            self.settings = {}
            
        logger.debug('self.settings: {}'.format(pformat(self.settings)))

        NetCDFConverter.__init__(self, 
                                 nc_out_path, 
                                 netcdf_format, 
                                 default_chunk_size=default_chunk_size, 
                                 default_variable_parameters=default_variable_parameters
                                 )
        
        self.aem_dat_path = aem_dat_path
        self.dfn_path = dfn_path
        
        self.field_definitions, self.crs = parse_dfn_file(dfn_path)
        
        # Read overriding field definition values from settings
        if self.settings.get('field_definitions'):
            for field_definition in self.field_definitions:
                overriding_field_definition = self.settings['field_definitions'].get(field_definition['short_name'])
                if overriding_field_definition:
                    field_definition.update(overriding_field_definition)
        
        logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
        
        #=======================================================================
        # # Determine index of 'nlayers' field definition. Field definitions after this will be 2D
        # self.max_dimension_field_index = self.field_definitions.index([field_def 
        #                                                        for field_def in self.field_definitions 
        #                                                        if field_def['short_name'] == 'layers'][0]
        #                                                        )
        # logger.debug('self.max_dimension_field_index: {}'.format(self.max_dimension_field_index))
        #=======================================================================
        
        self.dimension_field_definitions = {}
        for field_definition_index in range(len(self.field_definitions)):
            if self.field_definitions[field_definition_index]['short_name'] in self.settings['dimension_fields']:
                dimension_field_definition = dict(self.field_definitions[field_definition_index])
                
                self.dimension_field_definitions[field_definition_index] = dimension_field_definition
          
        logger.debug('self.dimension_field_definitions: {}'.format(pformat(self.dimension_field_definitions)))

        self.min_dimension_field_index = min(self.dimension_field_definitions.keys())
        self.max_dimension_field_index = max(self.dimension_field_definitions.keys())
        self.fields_per_dimension = len(self.field_definitions) - self.max_dimension_field_index - 1
        logger.debug('self.fields_per_dimension: {}'.format(self.fields_per_dimension))
        
        logger.info('Reading data file {}'.format(aem_dat_path))
        aem_dat_file = open(aem_dat_path, 'r')
         
        # Convert entire numeric file into a single array of float32 datatype
        # This is done to permit efficient column-wise data access but it will fail for any
        # invalid floating point values in data file
        #TODO: Use field widths from definition file
        row_list = []
        self.field_count = 0
        
        for line in aem_dat_file:
            line = line.strip()
            if not line: # Skip empty lines
                continue
            row = re.sub('\s+', ',', line).split(',')
            #logger.debug('row: {}'.format(row))
            
            if self.field_count:
                assert self.field_count == len(row), 'Inconsistent field count. Expected {}, found {}.'.format(self.field_count, 
                                                                                                               len(row)
                                                                                                               )
            else: # First row
                self.field_count = len(row)
            
            row_list.append(row)
            
        aem_dat_file.close()
         
        #logger.debug('row_list: {}'.format(row_list))
        self.raw_data_array = np.array(row_list, dtype='float32') # Convert list of lists to numpy array
        logger.info('{} points found'.format(self.raw_data_array.shape[0]))
        
        # Only check first row for secondary dimension sizes
        for dimension_field_index in self.dimension_field_definitions.keys():
            dimension_field_definition = self.dimension_field_definitions[dimension_field_index]
            dimension_name = dimension_field_definition['short_name']
            dimension_size = int(self.raw_data_array[0, dimension_field_index])
            dimension_field_definition['dimension_size'] = dimension_size 
            logger.info('secondary dimension "{}" is of size {}'.format(dimension_name, dimension_size))
        
        # Check layer_count for consistency across rows
        #TODO: Implement this check for multiple dimensions
        #assert not np.any(self.raw_data_array[:,self.max_dimension_field_index] - self.layer_count), 'Inconsistent layer count(s) found in column {}'.format(self.max_dimension_field_index + 1)
        
        # Check field count
        #TODO: Implement this check for multiple dimensions
        #expected_field_count = self.max_dimension_field_index + self.layer_count * self.fields_per_dimension + 1
        #assert self.field_count == expected_field_count, 'Invalid field count. Expected {}, found {}'.format(expected_field_count,
        #                                                                                                     self.field_count)
        
        # Process lines as a special case
        try:
            line_data = self.get_1d_data('line')
            
            self.lines, self.line_start_indices, self.line_point_counts = np.unique(line_data, 
                                                                                    return_index=True, 
                                                                                    return_inverse=False, 
                                                                                    return_counts=True
                                                                                    )
            logger.info('{} lines found'.format(len(self.lines)))
            #logger.debug('self.lines: {}, self.line_start_indices: {}, self.line_point_counts: {}'.format(self.lines, self.line_start_indices, self.line_point_counts))
        except:
            logger.info('No lines found')
            self.lines = None
            self.line_start_indices = None
            self.line_point_counts = None       
            
               
    def get_1d_data(self, short_name):
        '''
        Helper function to return 1D array corresponding to short_name from self.raw_data_array
        '''
        try:
            field_definition_index = self.field_definitions.index([field_def 
                                                        for field_def in self.field_definitions 
                                                        if field_def['short_name'] == short_name][0]
                                                        )
            return self.raw_data_array[:,field_definition_index]
        except:
            return None        
        
    def get_2d_data(self, short_name, dimension_name):
        '''
        Helper function to return 2D array corresponding to short_name from self.raw_data_array
        '''
        if True:#try:
            logger.debug('short_name: {}, dimension_name: {}'.format(short_name, dimension_name))
            field_definition_index = self.field_definitions.index([field_def 
                                                        for field_def in self.field_definitions 
                                                        if field_def['short_name'] == short_name][0]
                                                        )
            
            assert field_definition_index > self.max_dimension_field_index, 'Invalid 2D field name {}'.format(short_name)
            
            dimension_offset = field_definition_index - self.max_dimension_field_index - 1
            
            dimension_found = False
            field_count = 0
            field_offset = 0
            for dimension_field_index in sorted(self.dimension_field_definitions.keys()):
                logger.debug('dimension_field_index={}, self.dimension_field_definitions[dimension_field_index]={}'.format(dimension_field_index, self.dimension_field_definitions[dimension_field_index]))
                dimension_field_definition = self.dimension_field_definitions[dimension_field_index]
                dimension_size = dimension_field_definition['dimension_size']
                
                field_count += dimension_size
                
                if dimension_name == dimension_field_definition['short_name']:
                    dimension_found = True
                    break
                    
                field_offset += dimension_size
                    
            assert dimension_found, 'Dimension {} not found'.format(dimension_name)

            logger.debug('dimension_offset: {}, field_count: {}, field_offset: {}'.format(dimension_offset, field_count, field_offset))

            field_start_index = self.max_dimension_field_index + (dimension_offset * field_count) + field_offset + 1
            field_end_index = field_start_index + dimension_size
            
            logger.debug('field_start_index: {}, field_end_index: {}'.format(field_start_index, field_end_index))

            return self.raw_data_array[:,field_start_index:field_end_index]
        else:#except:
            return None        
        
    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs       
        '''
        metadata_dict = {'title': 'AEM .dat dataset read from {}'.format(os.path.basename(self.aem_dat_path)),
            'Conventions': "CF-1.6,ACDD-1.3",
            'featureType': "trajectory",
            'geospatial_east_min': np.min(self.get_1d_data('easting')),
            'geospatial_east_max': np.max(self.get_1d_data('easting')),
            'geospatial_east_units': "m",
            'geospatial_east_resolution': "point",
            'geospatial_north_min': np.min(self.get_1d_data('northing')),
            'geospatial_north_max': np.max(self.get_1d_data('northing')),
            'geospatial_north_units': "m",
            'geospatial_north_resolution': "point",
            'geospatial_vertical_min': np.min(self.get_1d_data('elevation')),
            'geospatial_vertical_max': np.max(self.get_1d_data('elevation')), # Should this be min(elevation-DOI)?
            'geospatial_vertical_units': "m",
            'geospatial_vertical_resolution': "point",
            'geospatial_vertical_positive': "up",
            'history': 'Converted from .dat file {} using definitions file {}'.format(self.aem_dat_path,
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
        dimensions = OrderedDict()
        
        # Example lat/lon dimensions
        dimensions['point'] = self.raw_data_array.shape[0]
        
        for dimension_field_index in sorted(self.dimension_field_definitions.keys()):
            dimension_field_definition = self.dimension_field_definitions[dimension_field_index]
            dimension_name = dimension_field_definition['short_name']
            dimension_size = dimension_field_definition['dimension_size']
            dimensions[dimension_name] = dimension_size 
            
        dimensions['line'] = len(self.lines)
              
        return dimensions
    
    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects       
        '''  
        def line_variable_generator():
            '''
            Helper generator to yield flight line variables - only called once
            '''
            assert self.lines is not None, 'flight Line parameters not defined'
            
            logger.info('\tWriting flight line indexing variables')
            
            logger.info('\t\tWriting flight line number')
            yield NetCDFVariable(short_name='lines', 
                                 data=self.lines, 
                                 dimensions=['line'], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'flight line number'}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting index of first point in flight line')
            yield NetCDFVariable(short_name='index_line', 
                                 data=self.line_start_indices, 
                                 dimensions=['line'], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'zero based index of the first sample in the line'}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
                        
            logger.info('\t\tWriting point count for flight line')
            yield NetCDFVariable(short_name='index_count', 
                                 data=self.line_point_counts, 
                                 dimensions=['line'], 
                                 fill_value=-1, 
                                 attributes={'long_name': 'number of samples in the line'}, 
                                 dtype='int32',
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
              
        # Create crs variable (points)
        yield self.build_crs_variable(self.crs)
        
        # Create 1D variables
        for field_1d_index in range(min(self.dimension_field_definitions.keys())):
            field_attributes = {}
            
            short_name = self.field_definitions[field_1d_index]['short_name']
            
            if short_name == 'line': # Special case for "line" variable
                for line_variable in line_variable_generator():
                    yield line_variable
                    
                continue
            
            long_name = self.field_definitions[field_1d_index].get('long_name')
            if long_name:
                field_attributes['long_name'] = long_name
                
            units = self.field_definitions[field_1d_index].get('units')
            if units:
                field_attributes['units'] = units
                
            fmt = self.field_definitions[field_1d_index].get('format')
            if fmt[0] =='I': # Integer field
                #TODO: see if we can reduce the size of the integer datatype
                dtype = 'int32' 
            else:
                dtype ='float32'
                
            fill_value=None
                
            logger.info('\tWriting 1D {} variable {}'.format(dtype, short_name))
                
            yield NetCDFVariable(short_name=short_name, 
                                 data=self.raw_data_array[:,field_1d_index], 
                                 dimensions=['point'], 
                                 fill_value=fill_value, 
                                 attributes=field_attributes, 
                                 dtype=dtype,
                                 chunk_size=self.default_chunk_size,
                                 variable_parameters=self.default_variable_parameters
                                 )
        
        #=======================================================================
        # # Create bad_data_mask array from depth of investigation
        # top_depth = self.get_2d_data('layer_top_depth')
        # depth_of_investigation = self.get_1d_data('depth_of_investigation')
        # bad_data_mask = top_depth > np.repeat(depth_of_investigation[:, np.newaxis], 
        #                                     top_depth.shape[1], 
        #                                     axis=1)
        # logger.debug('{} bad conductivity values found for masking'.format(np.count_nonzero(bad_data_mask)))
        #=======================================================================
        
        # Create 2D variables (points x <secondary dimension>)
        for field_2d_index in range(self.fields_per_dimension):
            field_definition_index = self.max_dimension_field_index + 1 + field_2d_index
            
            for dimension_field_index in sorted(self.dimension_field_definitions.keys()):
                dimension_field_definition = self.dimension_field_definitions[dimension_field_index]
                dimension_name = dimension_field_definition['short_name']
                dimension_long_name = dimension_field_definition['long_name']
                dimension_size = dimension_field_definition['dimension_size']

                field_attributes = {}
                
                short_name = self.field_definitions[field_definition_index]['short_name']
                
                logger.debug('field_2d_index: {}, field_definition_index: {}, dimension_field_index: {}'.format(field_2d_index, field_definition_index, dimension_field_index))
                
                long_name = self.field_definitions[field_definition_index].get('long_name')
                if long_name:
                    field_attributes['long_name'] = long_name + ' for ' + dimension_long_name
                     
                units = self.field_definitions[field_definition_index].get('units')
                if units:
                    field_attributes['units'] = units
                     
                fmt = self.field_definitions[field_definition_index].get('format')
                if fmt[0] =='I': # Integer field
                    #TODO: see if we can reduce the size of the integer datatype
                    dtype = 'int32' 
                else:
                    dtype='float32'
                logger.debug('fmt={}, dtype={}'.format(fmt, dtype))
                    
                # Convert resistivity to conductivity
                if short_name == 'resistivity':
                    short_name = 'conductivity'
                    field_attributes = {'long_name': 'Layer conductivity', 'units': 'S/m'}
                    data_array =  1.0 / self.get_2d_data('resistivity', dimension_name)
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
                    data_array =  self.get_2d_data('resistivity_uncertainty', dimension_name) / self.get_2d_data('resistivity', dimension_name)
                    fill_value = 0
                    #===========================================================
                    # data_array[bad_data_mask] = fill_value
                    #===========================================================
    
                    logger.debug('\nresistivity_uncertainty: {}'.format(self.get_2d_data('resistivity_uncertainty', dimension_name)))
                    logger.debug('\nresistivity: {}'.format(self.get_2d_data('resistivity', dimension_name)))
                    logger.debug('\nabsolute_resistivity_uncertainty: {}'.format(self.get_2d_data('resistivity', dimension_name) 
                                                                               * self.get_2d_data('resistivity_uncertainty', dimension_name)))
                    logger.debug('\nabsolute_conductivity_uncertainty: {}\n'.format(data_array))
                else:
                    # Don't mess with values
                    data_array = self.get_2d_data(short_name, dimension_name)
                    fill_value = None
                    
                # Create composite variable names reflecting dimensionality, e.g. conductivity_x_layers
                short_name += '_x_' + dimension_name
                if long_name:
                    long_name += ' for ' + dimension_long_name                
                        
                logger.info('\tWriting 2D {} variable {}'.format(dtype, short_name))
    
                yield NetCDFVariable(short_name=short_name, 
                                     data=data_array, 
                                     dimensions=['point', dimension_name], 
                                     fill_value=fill_value, 
                                     attributes=field_attributes, 
                                     dtype=dtype,
                                     chunk_size=self.default_chunk_size,
                                     variable_parameters=self.default_variable_parameters
                                     )
            
        return
    
def main():
    assert 3 <= len(sys.argv) <= 5, 'Invalid number of arguments.\n\
Usage: {} <dat_in_path> <dfn_in_path> [<nc_out_path>] [<settings_path>]'.format(sys.argv[0])
    dat_in_path = sys.argv[1] # 'C:\\Temp\\Groundwater Data\\ord_bonaparte_nbc_main_aquifer_clipped.dat'
    dfn_in_path = sys.argv[2] # 'C:\\Temp\\Groundwater Data\\nbc_20160421.dfn'

    if len(sys.argv) >= 4:
        nc_out_path = sys.argv[3] # 'C:\\Temp\\dat_test.nc'
    else:
        nc_out_path = os.path.splitext(dat_in_path)[0] + '.nc'
        
    if len(sys.argv) == 5:
        settings_path = sys.argv[4] # 'C:\\Temp\\dat_test.nc'
    else:
        settings_path = None
        
        
    d2n = AEMDAT2NetCDFConverter(nc_out_path, 
                                 dat_in_path, 
                                 dfn_in_path, 
                                 default_chunk_size=1024, 
                                 settings_path=settings_path)
    d2n.convert2netcdf()
    logger.info('Finished writing netCDF file {}'.format(nc_out_path))
    
    logger.debug('Global attributes:')
    logger.debug(pformat(d2n.nc_output_dataset.__dict__))
    logger.debug('Dimensions:')
    logger.debug(pformat(d2n.nc_output_dataset.dimensions))
    logger.debug('Variables:')
    for variable_name in d2n.nc_output_dataset.variables.keys():
        variable = d2n.nc_output_dataset.variables[variable_name]
        logger.debug(pformat(variable))
        print(variable[:])

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
