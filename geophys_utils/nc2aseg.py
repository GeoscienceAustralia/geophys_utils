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
Created on 8 Mar 2020

@author: Alex Ip <Alex.Ip@ga.gov.au>
'''

import argparse
import numpy as np
import re
import os
import sys
from datetime import datetime
from pprint import pformat
import tempfile
import logging
import locale
from math import log10
from collections import OrderedDict
from functools import reduce
import zipfile
import zipstream

from geophys_utils import get_spatial_ref_from_wkt
from geophys_utils import NetCDFPointUtils

locale.setlocale(locale.LC_ALL, '')  # Use '' for auto, or force e.g. to 'en_US.UTF-8'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

# Dynamically adjust integer field widths to fit all data values if True
ADJUST_INTEGER_FIELD_WIDTH = True

# Truncate ASEG-GDF2 field names to eight characters if True
TRUNCATE_VARIABLE_NAMES = False

# Set this to non zero to limit string field width in .dat file.
# WARNING - string truncation may corrupt data!
# N.B: Must be >= all numeric field widths defined in ASEG_GDF_FORMAT dict below (enforced by assertion)
MAX_FIELD_WIDTH = 0

# Maximum width of comment fields in .des file
MAX_COMMENT_WIDTH = 128

# Character encoding for .dfn, .dat & .des files
CHARACTER_ENCODING = 'utf-8'

# Default number of rows to read from netCDF before outputting a chunk of lines.
CACHE_CHUNK_ROWS = 32768
    
TEMP_DIR = tempfile.gettempdir()
#TEMP_DIR = 'C:\Temp'

# Set this to zero for no limit - only set a non-zero value for testing when debug = True
DEBUG_POINT_LIMIT = 100000

# List of regular expressions for variable names to exclude from output
EXCLUDE_NAME_REGEXES = ['crs', '.*_index$', 'ga_.*metadata', 'latitude.+', 'longitude.+', 'easting.+', 'northing.+']

# List of regular expressions for variable attributes to include in .dfn file
INCLUDE_VARIABLE_ATTRIBUTE_REGEXES = ['Intrepid.+']

# From Ross Brodie's email to Alex Ip, sent: Monday, 24 February 2020 4:27 PM
ASEG_GDF_FORMAT = {
    'float64': {
        'width': 18,
        'null': -9.9999999999e+32,
        'aseg_gdf_format': 'E18.10',
        'python_format': '{:>18.10e}',
        },
    'float32': {
        'width': 14,
        'null': -9.999999e+32,
        'aseg_gdf_format': 'E14.6',
        'python_format': '{:>14.6e}',
        },
    'int64': {
        'width': 21,
        'null': -9223372036854775808,
        'aseg_gdf_format': 'I21',
        'python_format': '{:>21d}',
        },
    'uint64': {
        'width': 21,
        'null': 18446744073709551616,
        'aseg_gdf_format': 'I21',
        'python_format': '{:>21d}',
        },
    'int32': {
        'width': 12,
        'null': -2147483647,
        'aseg_gdf_format': 'I12',
        'python_format': '{:>12d}',
        },
    'uint32': {
        'width': 12,
        'null': 4294967295,
        'aseg_gdf_format': 'I12',
        'python_format': '{:>12d}',
        },
    'int16': {
        'width': 7,
        'null': -32767,
        'aseg_gdf_format': 'I7',
        'python_format': '{:>7d}',
        },
    'uint16': {
        'width': 7,
        'null': 65535,
        'aseg_gdf_format': 'I7',
        'python_format': '{:>7d}',
        },
    'int8': {
        'width': 5,
        'null': -127,
        'aseg_gdf_format': 'I5',
        'python_format': '{:>5d}',
        },
    'uint8': {
        'width': 5,
        'null': 255,
        'aseg_gdf_format': 'I5',
        'python_format': '{:>5d}',
        },
    }

# Check to ensure that MAX_FIELD_WIDTH will not truncate numeric fields
assert not MAX_FIELD_WIDTH or all([format_specification['width'] <= MAX_FIELD_WIDTH for format_specification in ASEG_GDF_FORMAT.values()]), 'Invalid MAX_FIELD_WIDTH {}'.format(MAX_FIELD_WIDTH)

class RowValueCache(object):
    '''\
    Class to manage cache of row data from netCDF file
    '''
    def __init__(self, nc2aseggdf):
        '''
        Constructor
        '''
        self.nc2aseggdf = nc2aseggdf
        self.total_points = nc2aseggdf.total_points
        self.field_definitions = nc2aseggdf.field_definitions
        self.netcdf_dataset = nc2aseggdf.netcdf_dataset
        
        self.clear_cache()
       
        
    def clear_cache(self): 
        '''
        Clear cache
        '''
        self.index_range = 0
        self.cache = {}

    
    def read_points(self, start_index, end_index, point_mask=None):
        '''
        Function to read points from start_index to end_index
        '''
        self.index_range = end_index - start_index
        
        if point_mask is None: # No point_mask defined - take all points in range
            subset_mask = np.ones(shape=(self.index_range,), dtype='bool')
        else:
            subset_mask = point_mask[start_index:end_index]
            self.index_range = np.count_nonzero(subset_mask)
            
        # If no points to retrieve, don't read anything
        if not self.index_range:
            logger.debug('No points to retrieve - all masked out')            
            return

        # Build cache of data value slices keyed by field_name
        self.cache = {field_name: self.nc2aseggdf.get_data_values(field_name, slice(start_index, end_index))
                      for field_name in self.field_definitions.keys()
                      }
        
        #logger.debug('self.cache: {}'.format(pformat(self.cache)))
        
        
    def chunk_row_data_generator(self, clear_cache=True):
        '''
        Generator yielding chunks of all values from cache, expanding 2D variables to multiple columns
        '''
        if not self.index_range:
            logger.debug('Cache is empty - nothing to yield')
            return
        
        for index in range(self.index_range):
            row_value_list = []
            for field_name, field_definition in self.field_definitions.items():
                data = self.cache[field_name][index]
                
                # Convert array to string if required (OPeNDAP behaviour with string arrays?)
                if type(data) == np.ndarray and data.dtype == object:
                    data = str(data)

                if field_definition['columns'] == 1: # Element from 1D variable
                    row_value_list.append(data)
                elif field_definition['columns'] == 2: # Row from 2D variable
                    row_value_list += [element for element in data]
                else:
                    raise BaseException('Invalid dimensionality for variable {}'.format(field_definition['field_name']))
                            
            #logger.debug('row_value_list: {}'.format(row_value_list))
            yield row_value_list 
            
        if clear_cache:
            self.clear_cache() # Clear cache after outputting all lines                   
                                       
                                       
class NC2ASEGGDF2(object):
    
    def __init__(self,
                 netcdf_dataset, 
                 debug=False,
                 verbose=False,
                 ):
        '''
        '''
        def build_field_definitions():
            '''\
            Helper function to build self.field_definitions as an OrderedDict of field definitions keyed by ASEG-GDF2 field name
            '''
            self.field_definitions = OrderedDict()
            for variable_name, variable in self.netcdf_dataset.variables.items():
                
                # Check for any name exclusion matches
                if any([re.match(exclude_name_regex, variable_name, re.IGNORECASE)
                        for exclude_name_regex in EXCLUDE_NAME_REGEXES]):
                    logger.debug('Excluding variable {}'.format(variable_name))
                    continue
                
                if variable_name in self.field_definitions.keys(): # already processed
                    continue
                    
                if len(variable.dimensions) == 1 and variable.dimensions != ('point',): # Non-point indexed array variable, e.g.flight(line)
                    # Need to work backwards from variable to point indexed variable
                    try:
                        try:
                            index_variable = self.netcdf_dataset.variables[variable.index] # Try to use explicit index attribute
                        except AttributeError:
                            variable_dimension_name = variable.dimensions[0]
                            index_variable = self.netcdf_dataset.variables[variable_dimension_name + '_index']
                                    
                        assert index_variable.dimensions == ('point',), 'Invalid dimensions for variable {}: {}'.format(index_variable.name,
                                                                                                                        index_variable.dimensions)
                        
                        variables = [index_variable, variable]
                                                                                                                        
                        logger.debug('Found index variable {} for lookup variable {}'.format(index_variable.name, variable_name))
                    except:
                        logger.debug('Index variable not found for lookup variable {}'.format(variable_name))
                        continue # Can't find lookup variable - ignore this one   
                        
                elif (
                    len(variable.dimensions) 
                    and variable.dimensions[0] == 'point' 
                    and not (
                        variable.dimensions == ('point',) 
                        and (
                            variable_name.endswith('_index') 
                            or hasattr(variable, 'lookup')
                            )
                        )
                    ):
                    logger.debug('Found point-wise array data variable {}'.format(variable_name))
                    variables = [variable] # Not an index variable - just use primary variable values
                    
                elif not len(variable.dimensions) and variable_name != self.ncpu.crs_variable.name:
                    logger.debug('Found point-wise scalar data variable {}'.format(variable_name))
                    variables = [variable] # Scalar variable - broadcast out to all points
                    
                else:
                    logger.debug('Unable to deal with variable {} - ignoring'.format(variable_name))
                    continue
                
                variable_attributes = dict(variables[-1].__dict__)
                #logger.debug('variable_attributes = {}'.format(pformat(variable_attributes)))
                
                dtype = variables[-1].dtype
                logger.debug('Variable is of dtype {}'.format(dtype))
                format_dict = dict(ASEG_GDF_FORMAT.get(str(dtype)) or {})
                
                if not format_dict: # Unrecognised format. Treat as string
                    width = max([len(str(element).strip()) for element in variables[-1][:]]) + 1
                    
                    if MAX_FIELD_WIDTH and width > MAX_FIELD_WIDTH:
                        logger.warning('WARNING: String variable "{}" data will be truncated from a width of {} to {}'.format(variable_name, width, MAX_FIELD_WIDTH))
                        width = MAX_FIELD_WIDTH
                        
                    format_dict = {
                        'width': width,
                        'null': '',
                        'aseg_gdf_format': 'A{}'.format(width),
                        'python_format': '{{:>{}s}}'.format(width),
                        }
                    
                try:
                    column_count = reduce(lambda x, y: x*y, variable.shape[1:]) # This will work for (2+)D, even if we only support 1D or 2D
                except: # Scalar or 1D
                    column_count = 1
                
                variable_definition = {
                    'variable_name': variable_name,
                    'variables': variables,
                    'attributes': variable_attributes,
                    'dtype': dtype,
                    'format': format_dict,
                    'columns': column_count
                    }
                
                if TRUNCATE_VARIABLE_NAMES:
                    # Sanitise field name, truncate to 8 characters and ensure uniqueness
                    field_name = re.sub('(\W|_)+', '', variable_name)[:8].upper()
                    field_name_count = 0
                    while field_name in [variable_definition.get('field_name') 
                                         for variable_definition in self.field_definitions.values()]:
                        field_name_count += 1
                        field_name = field_name[:-len(str(field_name_count))] + str(field_name_count)
                else:
                    field_name = re.sub('\W+', '_', variable_name) # Sanitisation shouldn't be necessary, but we'll do it anyway
                    
                variable_definition['field_name'] = field_name   
                 
                # Add definition to allow subsequent self.get_data_values(field_name) call
                self.field_definitions[field_name] = variable_definition
                
                if ADJUST_INTEGER_FIELD_WIDTH and 'int' in str(dtype): # Field is some kind of integer - adjust format for data
                    #logger.debug('\tChecking values to adjust integer field width for variable {}'.format(variable_name))
                    max_abs_value = np.nanmax(np.abs(self.get_data_values(field_name)))
                    min_value = np.nanmin(self.get_data_values(field_name))
                    #logger.debug('\tMaximum absolute value = {}, minimum value = {}'.format(max_abs_value, min_value))
                    
                    if max_abs_value > 0:
                        width = int(log10(max_abs_value)) + 2 # allow for leading space
                        if min_value < 0:
                            width += 1 # allow for "-"
                    else:
                        width = 2
                        
                    if width != format_dict['width']:
                        logger.debug('\tAdjusting integer field width from {} to {} for variable {}'.format(format_dict['width'], width, variable_name))
                        format_dict['width'] = width
                        format_dict['aseg_gdf_format'] = 'I{}'.format(width)
                        format_dict['python_format'] = '{{:>{}d}}'.format(width)
                    
            #logger.debug(self.field_definitions)
        
        # Start of __init__
        
        #TODO: Make this a property
        self.debug = debug
        log_level = logging.DEBUG if debug else logging.INFO
        logger.setLevel(level=log_level)
        
        if verbose:
            logger.debug('Enabling info level output')
            self.info_output = logger.info # Verbose
        else:
            self.info_output = logger.debug # Non-verbose
                    
        self.ncpu = NetCDFPointUtils(netcdf_dataset, debug=debug)
        self.netcdf_dataset = self.ncpu.netcdf_dataset
        self.netcdf_path = self.ncpu.netcdf_dataset.filepath()
        
        self.netcdf_dataset.set_auto_mask(False) # Turn auto-masking off to allow substitution of new null values
        
        assert 'point' in self.netcdf_dataset.dimensions.keys(), '"point" not found in dataset dimensions'
        self.info_output('Opened netCDF dataset {}'.format(self.netcdf_path))
        
        self.total_points = self.ncpu.point_count
        
        self.spatial_ref = get_spatial_ref_from_wkt(self.ncpu.wkt)
        
        # set self.field_definitions 
        build_field_definitions()
        
        # set reporting increment to nice number giving 100 - 199 progress reports
        self.line_report_increment = (10.0 ** int(log10(self.ncpu.point_count / 50))) / 2.0
        
    
    def get_data_values(self, field_name, point_slice=slice(None, None, None)):
        '''\
        Function to return data values as an array, expanding lookups and broadcasting scalars if necessary
        @param field_name: Variable name to query (key in self.field_definitions)
        @param point_slice: slice to apply to point (i.e. first) dimension
        @return data_array: Array of data values
        '''
        variables = self.field_definitions[field_name]['variables'] 
        #logger.debug('Field {} represents variable {}({})'.format(field_name, variables[-1].name, ','.join(variables[0].dimensions)))
        
        if len(variables) == 1: # Single variable => no lookup
            if len(variables[0].dimensions): # Array 
                assert variables[0].dimensions[0] == 'point', 'First array dimension must be "point"'
                data = variables[0][point_slice]
            else: # Scalar
                # Broadcast scalar to required array shape
                data = np.array([variables[0][:]] * (((point_slice.stop or self.total_points) - (point_slice.start or 0)) // (point_slice.step or 1)))
        elif len(variables) == 2: # Index & Lookup variables
            data = variables[1][:][variables[0][:][point_slice]] # Use first array to index second one    
        else:
            raise BaseException('Unable to resolve chained lookups (yet): {}'.format([variable.name for variable in variables]))
        
        # Substitute null_value for _FillValue if required
        null_value = self.field_definitions[field_name]['format']['null']
        if null_value is not None and hasattr(variables[-1], '_FillValue'):
            data[(data == (variables[-1]._FillValue))] = null_value
        
        return data
    
    
    def create_dfn_line(
        self, 
        rt, 
        name, 
        aseg_gdf_format, 
        definition=None, 
        defn=None, 
        st='RECD'
        ):
        '''
        Helper function to write line to .dfn file.
        self.defn is used to track the DEFN number, which can be reset using the optional defn parameter
        @param rt: value for "RT=<rt>" portion of DEFN line, e.g. '' or 'PROJ'  
        @param name: Name of DEFN 
        @param format_specifier_dict: format specifier dict, e.g. {'width': 5, 'null': 256, 'aseg_gdf_format': 'I5', 'python_format': '{:>5d}'}
        @param definition=None: Definition string
        @param defn=None: New value of DEFN number. Defaults to self.defn+1
        @param st: value for "RT=<rt>" portion of DEFN line. Default = 'RECD'
        
        @return line: output line 
        '''
        if defn is None:
            self.defn += 1 # Increment last DEFN value (initialised to 0 in constructor)
        else:
            self.defn = defn
            
        line = 'DEFN {defn} ST={st},RT={rt}; {name}'.format(defn=self.defn,
                                                            st=st,
                                                            rt=rt,
                                                            name=name,
                                                            )    
        
        if aseg_gdf_format:
            line +=  ': {aseg_gdf_format}'.format(aseg_gdf_format=aseg_gdf_format)
            
        if definition:
            line += ': ' + definition
            
        #logger.debug('dfn file line: {}'.format(line))
        return line
        
                
    def create_dfn_file(self, dfn_out_path, zipstream_zipfile=None):
        '''
        Helper function to output .dfn file
        '''
        if zipstream_zipfile:
            dfn_basename = os.path.basename(dfn_out_path)
            
            zipstream_zipfile.write_iter(dfn_basename,
                                         self.encoded_dfn_line_generator(encoding=CHARACTER_ENCODING))
                    
        else:
            # Create, write and close .dfn file
            with open(dfn_out_path, 'w') as dfn_file:
                for dfn_line in self.dfn_line_generator():
                    dfn_file.write(dfn_line)
                dfn_file.close()
            self.info_output('Finished writing .dfn file {}'.format(self.dfn_out_path))
                
                    
    def encoded_dfn_line_generator(self, encoding=CHARACTER_ENCODING):
        '''
        Helper generator to yield encoded bytestrings of all lines in .dfn file
        '''
        for line_string in self.dfn_line_generator():
            yield line_string.encode(encoding) 
            
    def dfn_line_generator(self):
        '''
        Helper generator to yield all lines in .dfn file
        '''
        def variable_defns_generator():
            """
            Helper function to write a DEFN line for each variable
            """
            self.defn = 0 # reset DEFN number
            #for variable_name, variable_attributes in self.field_definitions.items():
            
            for field_name, field_definition in self.field_definitions.items():
                optional_attribute_list = []
                
                units = field_definition['attributes'].get('units')
                if units:
                    optional_attribute_list.append('UNITS={units}'.format(units=units))

                #fill_value = field_definition['attributes'].get('_FillValue')
                null_value = field_definition['format'].get('null')
                if null_value is not None:
                    optional_attribute_list.append('NULL=' + field_definition['format']['python_format'].format(null_value).strip())                   

                long_name = field_definition['attributes'].get('long_name') or re.sub('(\W|_)+', ' ', field_definition['variable_name'])
                if long_name:
                    optional_attribute_list.append('NAME={long_name}'.format(long_name=long_name))
                    
                # Include any variable attributes which match regexes in INCLUDE_VARIABLE_ATTRIBUTE_REGEXES
                for attribute_name, attribute_value in field_definition['attributes'].items():
                    if any([re.match(variable_attribute_regex, attribute_name, re.IGNORECASE)
                            for variable_attribute_regex in INCLUDE_VARIABLE_ATTRIBUTE_REGEXES]):
                        optional_attribute_list.append('{}={}'.format(attribute_name,
                                                                      attribute_value))

                #===========================================================
                # # Check for additional ASEG-GDF attributes defined in settings
                # variable_attributes = field_definition.get('variable_attributes')
                # if variable_attributes:
                #     for aseg_gdf_attribute, netcdf_attribute in self.settings['attributes'].items():
                #         attribute_value = variable_attributes.get(netcdf_attribute)
                #         if attribute_value is not None:
                #             optional_attribute_list.append('{aseg_gdf_attribute}={attribute_value}'.format(aseg_gdf_attribute=aseg_gdf_attribute,
                #                                                                                        attribute_value=attribute_value
                #                                                                                        ))
                #===========================================================
                    
                if optional_attribute_list:
                    definition = ', '.join(optional_attribute_list)
                else:
                    definition = None

                aseg_gdf_format = field_definition['format']['aseg_gdf_format']
                if field_definition['columns'] > 1: # Need to pre-pend number of columns to format string
                    aseg_gdf_format = '{}{}'.format(field_definition['columns'], aseg_gdf_format)
                
                yield self.create_dfn_line(rt='',
                                           name=field_name,
                                           aseg_gdf_format=aseg_gdf_format,
                                           definition=definition,
                                           )
                
                
            # Write 'END DEFN'
            yield self.create_dfn_line(rt='',
                                       name='END DEFN',
                                       aseg_gdf_format=None
                                       )
        
            
        def proj_defns_generator():
            """
            Helper function to write PROJ lines
From standard:
DEFN 1 ST=RECD,RT=PROJ; RT: A4
DEFN 2 ST=RECD,RT=PROJ; COORDSYS: A40: NAME=projection name, POSC projection name
DEFN 3 ST=RECD,RT=PROJ; DATUM: A40: NAME=datum name, EPSG compliant ellipsoid name
DEFN 4 ST=RECD,RT=PROJ; MAJ_AXIS: D12.1: UNIT=m, NAME=major_axis, Major axis in units
relevant to the ellipsoid definition
DEFN 5 ST=RECD,RT=PROJ; INVFLATT: D14.9: NAME=inverse flattening, 1/f inverse of flattening
DEFN 6 ST=RECD,RT=PROJ; PRIMEMER: F10.1: UNIT=deg, NAME=prime_meridian, Location of prime
meridian relative to Greenwich
DEFN 7 ST=RECD,RT=PROJ; PROJMETH: A30: NAME=projection_method, eg. Transverse Mercator,
Lambert etc
DEFN 8 ST=RECD,RT=PROJ; PARAM1: D14.0: NAME=Proj_par1, 1st projecton paramater See Table 1
DEFN 9 ST=RECD,RT=PROJ; PARAM2: D14.0: NAME=Proj_par2, 2nd projection parameter
DEFN 10 ST=RECD,RT=PROJ; PARAM3: D14.0: NAME=Proj_par3, 3rd projection parameter
DEFN 11 ST=RECD,RT=PROJ; PARAM4: D14.0: NAME=Proj_par4, 4th projection parameter
DEFN 12 ST=RECD,RT=PROJ; PARAM5: D14.0: NAME=Proj_par5, 5th projection parameter
DEFN 13 ST=RECD,RT=PROJ; PARAM6: D14.0: NAME=Proj_par6, 6th projection parameter
DEFN 14 ST=RECD,RT=PROJ; PARAM7: D14.0: NAME=Proj_par7, 7th projection parameter
DEFN 15 ST=RECD,RT=PROJ; END DEFN  

From sample file:
DEFN 1 ST=RECD,RT=PROJ; RT:A4
DEFN 2 ST=RECD,RT=PROJ; PROJNAME:A30: COMMENT=GDA94 / MGA zone 54
DEFN 3 ST=RECD,RT=PROJ; ELLPSNAM:A30: COMMENT=GRS 1980
DEFN 4 ST=RECD,RT=PROJ; MAJ_AXIS: D12.1: UNIT=m, COMMENT=6378137.000000
DEFN 5 ST=RECD,RT=PROJ; ECCENT: D12.9: COMMENT=298.257222
DEFN 6 ST=RECD,RT=PROJ; PRIMEMER: F10.1: UNIT=deg, COMMENT=0.000000
DEFN 7 ST=RECD,RT=PROJ; PROJMETH: A30: COMMENT=Transverse Mercator
DEFN 8 ST=RECD,RT=PROJ; PARAM1: D14.0: COMMENT=      0.000000
DEFN 9 ST=RECD,RT=PROJ; PARAM2: D14.0: COMMENT=    141.000000
DEFN 10 ST=RECD,RT=PROJ; PARAM3: D14.0: COMMENT=      0.999600
DEFN 11 ST=RECD,RT=PROJ; PARAM4: D14.0: COMMENT= 500000.000000
DEFN 12 ST=RECD,RT=PROJ; PARAM5: D14.0: COMMENT=10000000.00000
DEFN 13 ST=RECD,RT=PROJ; PARAM6: D14.0:
DEFN 14 ST=RECD,RT=PROJ; PARAM7: D14.0:
DEFN 15 ST=RECD,RT=PROJ; END DEFN
PROJGDA94 / MGA zone 54 GRS 1980  6378137.0000  298.257222  0.000000  Transverse Mercator  0.000000  141.000000  0.999600 500000.000000 10000000.00000
            """
            geogcs = self.spatial_ref.GetAttrValue('geogcs') # e.g. 'GDA94'
            projcs = self.spatial_ref.GetAttrValue('projcs') # e.g. 'UTM Zone 54, Southern Hemisphere'
            ellipse_name = self.spatial_ref.GetAttrValue('spheroid', 0)
            major_axis = float(self.spatial_ref.GetAttrValue('spheroid', 1))
            prime_meridian = float(self.spatial_ref.GetAttrValue('primem', 1))
            inverse_flattening = float(self.spatial_ref.GetInvFlattening())
            #eccentricity = self.spatial_ref.GetAttrValue('spheroid', 2) # Non-standard definition same as inverse_flattening?
            
            if self.spatial_ref.IsProjected():
                if projcs.startswith(geogcs):
                    projection_name = projcs
                else:
                    projection_name = geogcs + ' / ' + re.sub('[\:\,\=]+', '', projcs) # e.g. 'GDA94 / UTM Zone 54, Southern Hemisphere'
                projection_method = self.spatial_ref.GetAttrValue('projection').replace('_', ' ')
                projection_parameters = [(key, float(value))
                                          for key, value in re.findall('PARAMETER\["(.+)",(\d+\.?\d*)\]', self.spatial_ref.ExportToPrettyWkt())
                                          ]
            else: # Unprojected CRS
                projection_name = geogcs
                projection_method = None
                projection_parameters = None

            self.defn = 0  # reset DEFN number
            
            # write 'DEFN 1 ST=RECD,RT=PROJ; RT:A4'
            yield self.create_dfn_line(rt='PROJ',
                                       name='RT',
                                       aseg_gdf_format='A4'
                                       )
            
            yield self.create_dfn_line(rt='PROJ',
                                       name='COORDSYS',
                                       aseg_gdf_format='A40',
                                       definition='NAME={projection_name}, Projection name'.format(projection_name=projection_name)
                                       )

            yield self.create_dfn_line(rt='PROJ',
                                       name='DATUM',
                                       aseg_gdf_format='A40',
                                       definition='NAME={ellipse_name}, Ellipsoid name'.format(ellipse_name=ellipse_name)
                                       )
            
            yield self.create_dfn_line(rt='PROJ',
                                       name='MAJ_AXIS',
                                       aseg_gdf_format='D12.1',
                                       definition='UNIT={unit}, NAME={major_axis}, Major axis'.format(unit='m', major_axis=major_axis)
                                       )


            yield self.create_dfn_line(rt='PROJ',
                                       name='INVFLATT',
                                       aseg_gdf_format='D14.9',
                                       definition='NAME={inverse_flattening}, 1/f inverse of flattening'.format(inverse_flattening=inverse_flattening)
                                       )

            yield self.create_dfn_line(rt='PROJ',
                                       name='PRIMEMER',
                                       aseg_gdf_format='F10.1',
                                       definition='UNIT={unit}, NAME={prime_meridian}, Location of prime meridian'.format(unit='degree', prime_meridian=prime_meridian)
                                       )

#===============================================================================
#                 # Non-standard definitions
#                 yield self.create_dfn_line(rt='PROJ',
#                                            name='ELLPSNAM',
#                                            aseg_gdf_format='A30',
#                                            definition='NAME={ellipse_name}, Non-standard definition for ellipse name'.format(ellipse_name=ellipse_name)
#                                            )
# 
#                 yield self.create_dfn_line(rt='PROJ',
#                                            name='PROJNAME',
#                                            aseg_gdf_format='A40',
#                                            definition='NAME={projection_name}, Non-standard definition for projection name'.format(projection_name=projection_name)
#                                            )
# 
#                 yield self.create_dfn_line(rt='PROJ',
#                                            name='ECCENT',
#                                            aseg_gdf_format='D12.9',
#                                            definition='NAME={eccentricity}, Non-standard definition for ellipsoidal eccentricity'.format(eccentricity=eccentricity)
#                                            )
#===============================================================================
            
            if projection_method:                    
                yield self.create_dfn_line(rt='PROJ',
                                           name='PROJMETH',
                                           aseg_gdf_format='A30',
                                           definition='NAME={projection_method}, projection method'.format(projection_method=projection_method)
                                           )

                # Write all projection parameters starting from DEFN 8
                param_no = 0
                for param_name, param_value in projection_parameters:
                    param_no += 1  
                    yield self.create_dfn_line(rt='PROJ',
                                               name='PARAM{param_no}'.format(param_no=param_no),
                                               aseg_gdf_format='D14.0', #TODO: Investigate whether this is OK - it looks dodgy to me
                                               definition='NAME={param_value}, {param_name}'.format(param_value=param_value, param_name=param_name)
                                               )
            # Write 'END DEFN'
            yield self.create_dfn_line(rt='PROJ',
                                       name='END DEFN',
                                       aseg_gdf_format=''
                                       )
            
            #TODO: Write fixed length PROJ line at end of file
              
            return # End of function proj_defns_generator
            

        yield 'DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A{}\n'.format(MAX_COMMENT_WIDTH) # TODO: Check this first line 
        
        for defn_line in variable_defns_generator():
            yield defn_line + '\n'
        
        for proj_line in proj_defns_generator():
            yield proj_line + '\n'
    
    
    def create_dat_file(self, dat_out_path, cache_chunk_rows=None, point_mask=None, zipstream_zipfile=None):
        '''
        Helper function to output .dat file
        '''
        def chunk_buffer_generator(row_value_cache, python_format_list, cache_chunk_rows, point_mask=None, encoding=None):
            '''
            Generator to yield all line strings across all point variables for specified row range
            '''                    
            def chunk_line_generator(row_value_cache, python_format_list, start_index, end_index, point_mask=None):
                '''
                Helper Generator to yield line strings for specified rows across all point variables
                '''
                logger.debug('Reading rows {:n} - {:n}'.format(start_index+1, end_index))
                row_value_cache.read_points(start_index, end_index, point_mask=point_mask)
                
                logger.debug('Preparing ASEG-GDF lines for rows {:n} - {:n}'.format(start_index+1, end_index))
                for row_value_list in row_value_cache.chunk_row_data_generator():
                    #logger.debug('row_value_list: {}'.format(row_value_list))
                    # Turn list of values into a string using python_formats
                    # Truncate fields to maximum width with leading space - only string fields should be affected
                    yield ''.join([' ' + python_format_list[value_index].format(row_value_list[value_index])[1-MAX_FIELD_WIDTH::]
                                   for value_index in range(len(python_format_list))]) #.lstrip() # lstrip if we want to discard leading spaces from line
                    
            # Process all chunks
            point_count = 0
            for chunk_index in range(self.total_points // cache_chunk_rows + 1):
                chunk_line_list = []
                for line in chunk_line_generator(row_value_cache, python_format_list, start_index=chunk_index*cache_chunk_rows,
                                                 end_index=min((chunk_index+1)*cache_chunk_rows,
                                                               self.total_points
                                                               ),
                                                 point_mask=point_mask
                                                 ):
                    point_count += 1
                    
                    if not (point_count % self.line_report_increment):
                        self.info_output('{:n} / {:n} rows written'.format(point_count, self.total_points))
                    
                    #logger.debug('line: "{}"'.format(line))
                    chunk_line_list.append(line)
                
                    if self.debug and DEBUG_POINT_LIMIT and (point_count >= DEBUG_POINT_LIMIT): # Don't process more lines
                        break                    
                    
                chunk_buffer_string = '\n'.join(chunk_line_list) + '\n' # Yield a chunk of lines    
                
                if encoding:
                    yield(chunk_buffer_string.encode(encoding))
                else:
                    yield(chunk_buffer_string)
                        
                if self.debug and DEBUG_POINT_LIMIT and (point_count >= DEBUG_POINT_LIMIT): # Don't process more chunks
                    logger.warning('WARNING: Output limited to {:n} points in debug mode'.format(DEBUG_POINT_LIMIT))
                    break
                            
            self.info_output('A total of {:n} rows were output'.format(point_count))
            
            
        # Start of create_dat_file function
        cache_chunk_rows = cache_chunk_rows or CACHE_CHUNK_ROWS
        
        # Start of chunk_buffer_generator
        row_value_cache = RowValueCache(self) # Create cache for multiple chunks of data
        
        python_format_list = []
        for field_definition in self.field_definitions.values():
            for _column_index in range(field_definition['columns']):
                python_format_list.append(field_definition['format']['python_format'])
        #logger.debug('python_format_list: {}'.format(python_format_list)) 
                   
        if zipstream_zipfile:
            # Write to zip file
            dat_basename = os.path.basename(dat_out_path)
            zipstream_zipfile.write_iter(dat_basename, 
                         chunk_buffer_generator(row_value_cache, python_format_list, cache_chunk_rows, point_mask, encoding=CHARACTER_ENCODING))
        else: # No zip
            # Create, write and close .dat file
            dat_out_file = open(dat_out_path, mode='w')
            logger.debug('Writing lines to .dat file {}'.format(self.dat_out_path))
            for chunk_buffer in chunk_buffer_generator(row_value_cache, python_format_list, cache_chunk_rows, point_mask):
                logger.debug('Writing chunk_buffer to file')
                dat_out_file.write(chunk_buffer + '\n')
            dat_out_file.close()
            self.info_output('Finished writing .dat file {}'.format(dat_out_path))

    
    def create_des_file(self, des_out_path, zipstream_zipfile=None):
        '''
        Helper function to output .des file
        '''
        def des_line_generator(encoding=None):
            '''
            Helper Generator to yield line strings for .des file
            '''
            # Ignore netCDF system attributes
            global_attributes_dict = {key: str(value).strip() 
                                      for key, value in self.netcdf_dataset.__dict__.items()
                                      if not key.startswith('_')
                                      }
            
            # Determine maximum key length for fixed field width
            max_key_length = max([len(key) for key in global_attributes_dict.keys()])
                
            logger.debug('global_attributes_dict = {}'.format(pformat(global_attributes_dict)))   
            
            for key, value in global_attributes_dict.items():
                key_string = (' {{:<{}s}} : '.format(max_key_length)).format(key) # Include leading space
                
                for value_line in value.split('\n'):
                    # Split long values into multiple lines. Need to be careful with leading & trailing spaces when reassembling
                    while value_line: 
                        comment_line = 'COMM{}{}'.format(key_string, value_line[:MAX_COMMENT_WIDTH-len(key_string)]) +'\n'
                        
                        if encoding:
                            yield comment_line.encode(encoding)
                        else:
                            yield comment_line
                            
                        value_line = value_line[MAX_COMMENT_WIDTH-len(key_string):]
                    
        if zipstream_zipfile:
            # Write to zip file
            des_basename = os.path.basename(des_out_path)
            zipstream_zipfile.write_iter(des_basename, 
                                         des_line_generator(encoding=CHARACTER_ENCODING))
        else: # No zip
            # Create, write and close .dat file
            des_out_file = open(des_out_path, mode='w')
            logger.debug('Writing lines to .des file {}'.format(self.dat_out_path))
            for des_line in des_line_generator():
                logger.debug('Writing "{}" to .des file'.format(des_line))
                des_out_file.write(des_line)
            des_out_file.close()
            self.info_output('Finished writing .des file {}'.format(des_out_path))

    
    def convert2aseg_gdf(self, 
                         dat_out_path=None,
                         zip_out_path=None,
                         stride=1,
                         point_mask=None): 
        '''
        Function to convert netCDF file to ASEG-GDF
        '''
        self.dat_out_path = dat_out_path or os.path.splitext(self.netcdf_dataset.filepath())[0] + '.dat'
        self.dfn_out_path = os.path.splitext(dat_out_path)[0] + '.dfn'
        self.des_out_path = os.path.splitext(dat_out_path)[0] + '.des'
        
        if zip_out_path:
            zip_out_path = zip_out_path or os.path.splitext(dat_out_path)[0] + '_ASEG_GDF2.zip'
            zipstream_zipfile = zipstream.ZipFile(compression=zipfile.ZIP_DEFLATED)     
            zipstream_zipfile.comment = ('ASEG-GDF2 files generated at {} from {}'.format(datetime.now().isoformat(),
                                                                                          os.path.basename(self.netcdf_path))
                                                                                          ).encode(CHARACTER_ENCODING)  
            try:
                os.remove(zip_out_path)
            except:
                pass
            
        else:
            zipstream_zipfile = None

        self.create_dfn_file(self.dfn_out_path, zipstream_zipfile=zipstream_zipfile)

        self.create_dat_file(self.dat_out_path, zipstream_zipfile=zipstream_zipfile)
        
        self.create_des_file(self.des_out_path, zipstream_zipfile=zipstream_zipfile)
        
        if zipstream_zipfile:
            logger.debug('Opening zip file {}'.format(zip_out_path))
            zip_out_file = open(zip_out_path, 'wb')
            
            self.info_output('Writing zip file {}'.format(zip_out_path))
            for data in zipstream_zipfile:
                zip_out_file.write(data)
                
            self.info_output('Closing zip file {}'.format(zip_out_path))
            zipstream_zipfile.close()



def main():
    '''
    Main function
    '''
    def get_args():
        """
        Handles all the arguments that are passed into the script

        :return: Returns a parsed version of the arguments.
        """
        parser = argparse.ArgumentParser(description='Convert netCDF file to ASEG-GDF2')

        parser.add_argument("-r", "--crs",
                            help="Coordinate Reference System string (e.g. GDA94, EPSG:4283) for output",
                            type=str,
                            dest="crs")
        
        parser.add_argument('-z', '--zip', action='store_const', const=True, default=False,
                            help='Zip directly to an archive file. Default is no zip')
        
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='output debug information. Default is no debug info')
        
        parser.add_argument('-v', '--verbose', action='store_const', const=True, default=False,
                            help='output verbosity. Default is non-verbose')
        
        parser.add_argument('positional_args', 
                            nargs=argparse.REMAINDER,
                            help='<nc_in_path> [<dat_out_path>]')

        return parser.parse_args()
    
    args = get_args()
    
    # Setup Logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logger.setLevel(level=log_level)

    assert 1 <= len(args.positional_args) <= 2, 'Invalid number of positional arguments.\n\
Usage: python {} <options> <nc_in_path> [<dat_out_path>] [<zip_out_path>]'.format(os.path.basename(sys.argv[0]))

    nc_in_path = args.positional_args[0]

    if len(args.positional_args) == 2:
        dat_out_path = args.positional_args[1]
    else:
        dat_out_path = os.path.splitext(nc_in_path)[0] + '.dat'
        
    if args.zip:
        if len(args.positional_args) == 3:
            zip_out_path = args.positional_args[2]
        else:
            zip_out_path = os.path.splitext(nc_in_path)[0] + '.zip'
    else:
        zip_out_path = None
        
    logger.debug('args: {}'.format(args.__dict__))

    nc2aseggdf2 = NC2ASEGGDF2(nc_in_path, debug=args.debug, verbose=args.verbose)
    
    nc2aseggdf2.convert2aseg_gdf(dat_out_path, zip_out_path)

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