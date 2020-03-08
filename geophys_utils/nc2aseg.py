'''
Created on 8 Mar 2020

@author: alex
'''
import argparse
import numpy as np
import re
import os
import sys
from datetime import datetime
from pprint import pformat, pprint
import yaml
import tempfile
import netCDF4
import logging
import locale
from math import log10
from collections import OrderedDict

from geophys_utils.netcdf_converter.aseg_gdf_utils import variable2aseg_gdf_format
from geophys_utils import get_spatial_ref_from_wkt
from geophys_utils import NetCDFPointUtils

locale.setlocale(locale.LC_ALL, '')  # Use '' for auto, or force e.g. to 'en_US.UTF-8'

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

# From Ross Brodie's email to Alex Ip, sent: Monday, 24 February 2020 4:27 PM
ASEG_GDF_FORMAT = {
    'float64': {
        'width': 18,
        'null': -9.9999999999e+00,
        'aseg_gdf_format': 'E18.10',
        },
    'float32': {
        'width': 14,
        'null': -9.999999e+00,
        'aseg_gdf_format': 'E14.6',
        },
    'int32': {
        'width': 12,
        'null': -2147483647,
        'aseg_gdf_format': 'I12',
        },
    'uint32': {
        'width': 12,
        'null': 4294967296,
        'aseg_gdf_format': 'I12',
        },
    'int16': {
        'width': 7,
        'null': -32767,
        'aseg_gdf_format': 'I7',
        },
    'uint16': {
        'width': 7,
        'null': 65536,
        'aseg_gdf_format': 'I7',
        },
    'int8': {
        'width': 5,
        'null': -127,
        'aseg_gdf_format': 'I5',
        },
    'uint8': {
        'width': 5,
        'null': 256,
        'aseg_gdf_format': 'I5',
        },
    }

# Default number of rows to read before outputting chunk of lines.
CACHE_CHUNK_ROWS = 8192
    
TEMP_DIR = tempfile.gettempdir()
#TEMP_DIR = 'C:\Temp'

# Set this to zero for no limit - only set a non-zero value for testing
POINT_LIMIT = 1000000

EXCLUDE_NAME_REGEXES = ['latitude.+', 'longitude.+', 'easting.+', 'northing.+']

class NC2ASEGGDF2(object):
    
    def __init__(self,
                 netcdf_dataset, 
                 dat_file_path=None, 
                 dfn_file_path=None, 
                 debug=False):
        '''
        '''
        def build_variable_definitions():
            variable_definitions = OrderedDict()
            for variable_name, variable in self.netcdf_dataset.variables.items():
                
                # Check for any name exclusion matches
                if any([re.match(exclude_name_regex, variable_name)
                        for exclude_name_regex in EXCLUDE_NAME_REGEXES]):
                    logger.debug('Excluding variable {}'.format(variable_name))
                    continue
                
                if variable_name in variable_definitions.keys():
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
                    
                elif not len(variable.dimensions) and variable_name != ncpu.crs_variable.name:
                    logger.debug('Found point-wise scalar data variable {}'.format(variable_name))
                    variables = [variable] # Scalar variable - broadcast out to all points
                    
                else:
                    logger.debug('Unable to deal with variable {} - ignoring'.format(variable_name))
                    continue
                
                variable_attributes = dict(variables[-1].__dict__)
                #logger.debug('variable_attributes = {}'.format(pformat(variable_attributes)))
                
                dtype = variables[-1].dtype
                variable_definitions[variable_name] = {
                    'variables': variables,
                    'attributes': variable_attributes,
                    'dtype': dtype,
                    'format': ASEG_GDF_FORMAT[str(dtype)] # Don't use "get" so we flush out missing formats
                    }
                    
            return variable_definitions
        
        # Start of __init__
        dat_file_path = dat_file_path or os.path.splitext(netcdf_dataset)[0] + '.dat'
        dfn_file_path = dfn_file_path or os.path.splitext(dat_file_path)[0] + '.dfn'
        ncpu = NetCDFPointUtils(netcdf_dataset, debug=debug)
        self.netcdf_dataset = ncpu.netcdf_dataset
        
        assert 'point' in self.netcdf_dataset.dimensions.keys(), '"point" not found in dataset dimensions'
        
        self.variable_definitions = build_variable_definitions()
        #logger.debug('variable_definitions = {}'.format(pformat(self.variable_definitions)))    
        
        for variable_name in self.variable_definitions.keys():
            print('{} = {}'.format(variable_name, self.get_data_values(variable_name, slice(None, None, 100000))))
        
    
    def get_data_values(self, variable_name, point_slice=slice(None, None, None)):
        '''\
        Function to return data values as an array
        @param variable_name: Variable name to query (key in self.variable_definitions)
        @param point_slice: slice to apply to point (i.e. first) dimension
        @return data_array: Array of data values
        '''
        variables = self.variable_definitions[variable_name]['variables'] 
        #logger.debug('{}({})'.format(variable_name, ','.join(variables[0].dimensions)))
        
        if len(variables) == 1: # Single variable => no lookup
            if len(variables[0].dimensions): # Array 
                assert variables[0].dimensions[0] == 'point', 'First array dimension must be "point"'
                return variables[0][:][point_slice]
            else: # Scalar
                # Broadcast scalar to required array shape
                return np.array([variables[0][:]] * (((point_slice.stop or len(variables[0])) - (point_slice.start or 0)) // (point_slice.step or 1)))
        elif len(variables) == 2: # Index & Lookup variables
            return variables[1][:][variables[0][:][point_slice]] # Use first array to index second one    
        else:
            raise BaseException('Unable to resolve chained lookups: {}'.format([variable.name for variable in variables]))
        
    
    def write_record2dfn_file(self, 
                              dfn_file, 
                              rt, 
                              name, 
                              aseg_gdf_format, 
                              definition=None, 
                              defn=None, 
                              st='RECD'):
        '''
        Helper function to write line to .dfn file.
        self.defn is used to track the DEFN number, which can be reset using the optional defn parameter
        @param dfn_file: output file for DEFN line 
        @param rt: value for "RT=<rt>" portion of DEFN line, e.g. '' or 'PROJ'  
        @param name: Name of DEFN 
        @param aseg_gdf_format: ASEG-GDF output format, e.g. 'I5', 'D12.1' or 'A30'
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
            
        dfn_file.write(line + '\n')
        return line
        
                
    def convert2aseg_gdf(self, 
                         dat_out_path=None,
                         dfn_out_path=None,
                         point_mask=None): 
        '''
        Function to convert netCDF file to ASEG-GDF
        '''
        def get_field_definitions():
            '''
            Helper function to get field definitions from netCDF file
            '''
            # Build field definitions
            self.field_definitions = []
            for variable_name, variable in self.nc_dataset.variables.items():
                # Skip any non-scalar, non-pointwise fields except index variables
                # These are probably lookup fields
                if (variable.shape 
                    and not (
                        ('point' in variable.dimensions) # Point-indexed variables
                        or (variable.dimensions == ('line',) and variable_name != 'line') # Line-indexed variables
                        )):
                    continue
                
                # Don't output CRS scalar variable - assume all other scalars are to be applied to every point
                if (not variable.shape
                    and (variable_name in ['crs', 'transverse_mercator']
                         or re.match('ga_.+_metadata', variable_name)
                         )
                    ):
                        continue
                
                # Resolve lookups - use lookup variable instead of index variable
                if variable_name.endswith('_index') | hasattr(variable, 'lookup'): # Index variable found
                    # Use lookup variable name for output
                    try:
                        short_name = variable.lookup # Use lookup attribute
                    except AttributeError:
                        short_name = re.sub('_index$', '', variable_name, re.IGNORECASE) # Use associated variable name
                    
                    try:
                        lookup_variable = self.nc_dataset.variables[short_name]
                    except:
                        logger.warning('Unable to access lookup variable "{}" for index variable "{}"'.format(short_name, variable_name))
                        continue
                    
                    variable_attributes = lookup_variable.__dict__
                    aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format = variable2aseg_gdf_format(lookup_variable)

                    
                else: # Non-lookup or line-indexed variable
                    short_name = variable_name
                    variable_attributes = variable.__dict__
                    aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format = variable2aseg_gdf_format(variable)
                    
                # Map variable name to standard ASEG-GDF field name if required
                short_name = self.settings['aseg_field_mapping'].get(short_name) or short_name
                    
                try:
                    read_chunk_size = int(variable.chunking()[0])
                except:
                    read_chunk_size = CACHE_CHUNK_ROWS # Use default chunking for reads
                    
                field_definition = {'variable_name': variable_name,
                                    'short_name': short_name,
                                    'dtype': dtype,
                                    'chunk_size': read_chunk_size,
                                    'columns': columns,
                                    'format': aseg_gdf_format,
                                    'width_specifier': width_specifier,
                                    'decimal_places': decimal_places,
                                    'python_format': python_format
                                    }
                
                fill_value = variable_attributes.get('_FillValue') # System attribute
                if fill_value is not None:
                    field_definition['fill_value'] = fill_value
                                
                long_name = variable_attributes.get('long_name')
                if long_name is not None:
                    field_definition['long_name'] = long_name
                                
                # Set variable attributes in field definition
                variable_attribute_dict = {attribute_name: variable_attributes.get(key.upper())
                    for key, attribute_name in self.settings['variable_attributes'].items()
                    if variable_attributes.get(key.upper()) is not None
                                           }
                        
                if variable_attribute_dict:
                    field_definition['variable_attributes'] = variable_attribute_dict  
                      
                self.field_definitions.append(field_definition)
                
            logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
    
            # Read overriding field definition values from settings
            if self.settings.get('field_definitions'):
                for field_definition in self.field_definitions:
                    overriding_field_definition = self.settings['field_definitions'].get(field_definition['short_name'])
                    if overriding_field_definition:
                        field_definition.update(overriding_field_definition)
                
                logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
        
        
        def write_dfn_file():
            '''
            Helper function to output .dfn file
            '''
            def write_defns(dfn_file):
                """
                Helper function to write multiple DEFN lines
                """
                self.defn = 0 # reset DEFN number
                for field_definition in self.field_definitions:
                    short_name = field_definition['short_name']
                    
                    optional_attribute_list = []
                    
                    units = field_definition.get('units')
                    if units:
                        optional_attribute_list.append('UNITS={units}'.format(units=units))
    
                    fill_value = field_definition.get('fill_value')
                    if fill_value is not None:
                        optional_attribute_list.append('NULL=' + field_definition['python_format'].format(fill_value).strip())                   
    
                    long_name = field_definition.get('long_name')
                    if long_name:
                        optional_attribute_list.append('NAME={long_name}'.format(long_name=long_name))
                        
                    # Check for additional ASEG-GDF attributes defined in settings
                    variable_attributes = field_definition.get('variable_attributes')
                    if variable_attributes:
                        for aseg_gdf_attribute, netcdf_attribute in self.settings['variable_attributes'].items():
                            attribute_value = variable_attributes.get(netcdf_attribute)
                            if attribute_value is not None:
                                optional_attribute_list.append('{aseg_gdf_attribute}={attribute_value}'.format(aseg_gdf_attribute=aseg_gdf_attribute,
                                                                                                           attribute_value=attribute_value
                                                                                                           ))
                        
                    if optional_attribute_list:
                        definition = ' , '.join(optional_attribute_list)
                    else:
                        definition = None
    
                    self.write_record2dfn_file(dfn_file,
                                               rt='',
                                               name=short_name,
                                               aseg_gdf_format=field_definition['format'],
                                               definition=definition,
                                               )
                    
                    
                # Write 'END DEFN'
                self.write_record2dfn_file(dfn_file,
                                           rt='',
                                           name='END DEFN',
                                           aseg_gdf_format=''
                                           )
                    
                return # End of function write_defns
            
                
            def write_proj(dfn_file):
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
                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='RT',
                                           aseg_gdf_format='A4'
                                           )
                
                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='COORDSYS',
                                           aseg_gdf_format='A40',
                                           definition='NAME={projection_name}, Projection name'.format(projection_name=projection_name)
                                           )

                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='DATUM',
                                           aseg_gdf_format='A40',
                                           definition='NAME={ellipse_name}, Ellipsoid name'.format(ellipse_name=ellipse_name)
                                           )
                
                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='MAJ_AXIS',
                                           aseg_gdf_format='D12.1',
                                           definition='UNIT={unit}, NAME={major_axis}, Major axis'.format(unit='m', major_axis=major_axis)
                                           )


                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='INVFLATT',
                                           aseg_gdf_format='D14.9',
                                           definition='NAME={inverse_flattening}, 1/f inverse of flattening'.format(inverse_flattening=inverse_flattening)
                                           )

                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='PRIMEMER',
                                           aseg_gdf_format='F10.1',
                                           definition='UNIT={unit}, NAME={prime_meridian}, Location of prime meridian'.format(unit='degree', prime_meridian=prime_meridian)
                                           )

#===============================================================================
#                 # Non-standard definitions
#                 self.write_record2dfn_file(dfn_file,
#                                            rt='PROJ',
#                                            name='ELLPSNAM',
#                                            aseg_gdf_format='A30',
#                                            definition='NAME={ellipse_name}, Non-standard definition for ellipse name'.format(ellipse_name=ellipse_name)
#                                            )
# 
#                 self.write_record2dfn_file(dfn_file,
#                                            rt='PROJ',
#                                            name='PROJNAME',
#                                            aseg_gdf_format='A40',
#                                            definition='NAME={projection_name}, Non-standard definition for projection name'.format(projection_name=projection_name)
#                                            )
# 
#                 self.write_record2dfn_file(dfn_file,
#                                            rt='PROJ',
#                                            name='ECCENT',
#                                            aseg_gdf_format='D12.9',
#                                            definition='NAME={eccentricity}, Non-standard definition for ellipsoidal eccentricity'.format(eccentricity=eccentricity)
#                                            )
#===============================================================================
                
                if projection_method:                    
                    self.write_record2dfn_file(dfn_file,
                                               rt='PROJ',
                                               name='PROJMETH',
                                               aseg_gdf_format='A30',
                                               definition='NAME={projection_method}, projection method'.format(projection_method=projection_method)
                                               )

                    # Write all projection parameters starting from DEFN 8
                    param_no = 0
                    for param_name, param_value in projection_parameters:
                        param_no += 1  
                        self.write_record2dfn_file(dfn_file,
                                                   rt='PROJ',
                                                   name='PARAM{param_no}'.format(param_no=param_no),
                                                   aseg_gdf_format='D14.0', #TODO: Investigate whether this is OK - it looks dodgy to me
                                                   definition='NAME={param_value}, {param_name}'.format(param_value=param_value, param_name=param_name)
                                                   )
                # Write 'END DEFN'
                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='END DEFN',
                                           aseg_gdf_format=''
                                           )
                
                #TODO: Write fixed length PROJ line at end of file
                  
                return # End of function write_proj
                

            # Create, write and close .dat file
            dfn_file = open(self.dfn_out_path, 'w')
            dfn_file.write('DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76\n') # TODO: Check this first line 
            
            write_defns(dfn_file)
            
            write_proj(dfn_file)
                
            dfn_file.close()
            self.info_output('Finished writing .dfn file {}'.format(self.dfn_out_path))
        
        
        def write_dat_file(cache_chunk_rows=None, point_mask=None):
            '''
            Helper function to output .dat file
            '''
            def chunk_buffer_generator(point_mask=None):
                '''
                Generator to yield all line strings across all point variables for specified row range
                '''                    
                def chunk_line_generator(start_index, end_index, point_mask=None):
                    '''
                    Helper Generator to yield line strings for specified rows across all point variables
                    '''
                    python_format_list = []
                    for field_definition in self.field_definitions:
                        for _column_index in range(field_definition['columns']):
                            python_format_list.append(field_definition['python_format'])
                    logger.debug('python_format_list: {}'.format(python_format_list)) 
                           
                    value_count = len(python_format_list)
        
                    logger.debug('Reading rows {:n}-{:n}'.format(start_index+1, end_index))
                    line_cache.read_points(start_index, end_index, point_mask=point_mask)
                    
                    logger.debug('Preparing ASEG-GDF lines for rows {:n}-{:n}'.format(start_index+1, end_index))
                    for value_list in line_cache.chunk_buffer_generator():
                        logger.debug('value_list: {}'.format(value_list))
                        # Turn list of values into a string using python_formats
                        yield ''.join([python_format_list[value_index].format(value_list[value_index])
                                       for value_index in range(value_count)]).lstrip()

                        
                # Start of chunk_buffer_generator
                line_cache = RowCache(self) # Create cache for multiple rows
                
                # Process all chunks
                point_count = 0
                for chunk_index in range(self.total_points // cache_chunk_rows + 1):
                    for line in chunk_line_generator(start_index=chunk_index*cache_chunk_rows,
                                                     end_index=min((chunk_index+1)*cache_chunk_rows,
                                                                   self.total_points
                                                                   ),
                                                     point_mask=point_mask
                                                     ):
                        point_count += 1
                        
                        if point_count == point_count // self.line_report_increment * self.line_report_increment:
                            self.info_output('{:n} / {:n} rows written'.format(point_count, self.total_points))
                        
                        logger.debug('line: {}'.format(line))
                        yield line
                    
                        if POINT_LIMIT and (point_count >= POINT_LIMIT):
                            break                    
                        
                    if POINT_LIMIT and (point_count >= POINT_LIMIT):
                        break
                
                self.info_output('{:n} rows output'.format(point_count))
                
                
            # Start of write_dat_file function
            cache_chunk_rows = cache_chunk_rows or CACHE_CHUNK_ROWS
            
            # Create, write and close .dat file
            dat_out_file = open(self.dat_out_path, mode='w')
            logger.debug('Writing lines to {}'.format(self.dat_out_path))
            for line in chunk_buffer_generator(point_mask):
                dat_out_file.write(line + '\n')
            dat_out_file.close()
            self.info_output('Finished writing .dat file {}'.format(self.dat_out_path))
                
        
        # Start of convert2aseg_gdf function        
        self.dat_out_path = dat_out_path or os.path.splitext(self.netcdf_in_path)[0] + '.dat'
        self.dfn_out_path = dfn_out_path or os.path.splitext(dat_out_path)[0] + '.dfn'
        
        get_field_definitions()
        
        write_dfn_file()

        write_dat_file(point_mask=point_mask)
        
        
    

    



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
Usage: python {} <options> <nc_in_path> [<dat_out_path>]'.format(os.path.basename(sys.argv[0]))

    nc_in_path = args.positional_args[0]

    if len(args.positional_args) == 2:
        dat_out_path = args.positional_args[1]
    else:
        dat_out_path = os.path.splitext(nc_in_path)[0] + '.dat'
        
    dfn_out_path = os.path.splitext(dat_out_path)[0] + '.dfn'
    
    logger.debug('args: {}'.format(args.__dict__))

    nc2aseggdf2 = NC2ASEGGDF2(nc_in_path, dat_out_path, dfn_out_path)


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