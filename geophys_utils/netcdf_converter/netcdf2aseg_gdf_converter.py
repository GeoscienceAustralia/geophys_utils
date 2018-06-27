'''
Created on 13 Jun. 2018

@author: u76345
'''

import argparse
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

from geophys_utils.netcdf_converter.aseg_gdf_utils import variable2aseg_gdf_format
from geophys_utils import get_spatial_ref_from_wkt


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # Logging level for this module

# Maximum size of in-memory cache array (in bytes)
MAX_MEMORY_BYTES = 1073741824 # 1GB
#MAX_MEMORY_BYTES = 8388608 # 8MB

# Default effective chunk size for un-chunked variables.
# Set to 0 for no chunking (i.e. complete read)
DEFAULT_READ_CHUNK_SIZE = 1
    
TEMP_DIR = tempfile.gettempdir()
#TEMP_DIR = 'C:\Temp'

# Set this to zero for no limit - only set a non-zero value for testing
POINT_LIMIT = 0

class NetCDF2ASEGGDFConverter(object):
    '''
    NetCDF2ASEGGDFConverter class definition to convert netCDF file to ASEG-GDF format
    '''
    def __init__(self,
                 netcdf_in_path,
                 dat_out_path=None,
                 dfn_out_path=None,
                 settings_path=None
                 ):
        '''
        Constructor for class NetCDF2ASEGGDFConverter
        '''
        assert os.path.isfile(netcdf_in_path), '{} is not a valid file'.format(netcdf_in_path)
        
        self.netcdf_in_path = netcdf_in_path
        self.dat_out_path = dat_out_path or os.path.splitext(netcdf_in_path)[0] + '.dat'
        self.dfn_out_path = dfn_out_path or os.path.splitext(dat_out_path)[0] + '.dfn'
        self.defn = 0 # DEFN number for last line written to .dfn file
        
        self.settings_path = settings_path or os.path.join(os.path.dirname(__file__), 
                                                           'aseg_gdf_settings.yml')
        
        try:
            self.settings = yaml.safe_load(open(self.settings_path))
        except:
            self.settings = {}
            
        logger.debug('self.settings: {}'.format(pformat(self.settings)))

        self.nc_dataset = netCDF4.Dataset(self.netcdf_in_path, 'r')
        assert 'point' in self.nc_dataset.dimensions.keys(), 'No point dimension defined in netCDF dataset'
        
        self.total_points = self.nc_dataset.dimensions['point'].size
        
        # Build field definitions
        self.field_definitions = []
        for variable_name, variable in self.nc_dataset.variables.items():
            
            # Skip scalar variables
            if not len(variable.dimensions):
                continue
            
            # Skip any non-pointwise, non-indexing fields
            # Note: Indexing variables will be expanded during .dat output
            if ('point' not in variable.dimensions) and variable_name not in self.settings['index_fields']:
                continue
            
            chunk_size = variable.chunking()[0]
            # If variable is not chunked and DEFAULT_READ_CHUNK_SIZE is defined
            if DEFAULT_READ_CHUNK_SIZE and (chunk_size == variable.shape[0]):
                chunk_size = min(DEFAULT_READ_CHUNK_SIZE, variable.shape[0]) # Use default chunking
                
            aseg_gdf_format, dtype, columns, integer_digits, fractional_digits, python_format = variable2aseg_gdf_format(variable)
            
            #TODO: Add extra field definition stuff like ASEG-GDF format specifier
            field_definition = {'short_name': variable_name,
                                'dtype': dtype,
                                'chunk_size': chunk_size,
                                'columns': columns,
                                'format': aseg_gdf_format,
                                'integer_digits': integer_digits,
                                'fractional_digits': fractional_digits,
                                'python_format': python_format
                                }
            
            self.field_definitions.append(field_definition)

        logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))

        # Read overriding field definition values from settings
        if self.settings.get('field_definitions'):
            print(self.settings['field_definitions'])
            for field_definition in self.field_definitions:
                print('short_name:', field_definition['short_name'])
                overriding_field_definition = self.settings['field_definitions'].get(field_definition['short_name'])
                if overriding_field_definition:
                    field_definition.update(overriding_field_definition)
            
            logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
    
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
        
                
    def convert2aseg_gdf(self): 
        '''
        Function to convert netCDF file to ASEG-GDF
        '''
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
                    variable_name = field_definition['short_name']
                    # Map variable name to standard ASEG-GDF field name if required
                    field_name = self.settings['aseg_field_mapping'].get(variable_name) or variable_name
                    
                    variable = self.nc_dataset.variables[variable_name]
                    logger.debug('variable_name: {}, field_name: {}, variable.dtype: {}'.format(variable_name, field_name, variable.dtype))
                    
                    variable_attributes = variable.__dict__
                    
                    optional_attribute_list = []
                    
                    units = variable_attributes.get('units')
                    if units:
                        optional_attribute_list.append('UNITS={units}'.format(units=units))
    
                    fill_value = variable_attributes.get('_FillValue')
                    if fill_value is not None:
                        optional_attribute_list.append('NULL=' + field_definition['python_format'].format(fill_value).strip())                   
    
                    # Check for additional ASEG-GDF attributes defined in settings
                    for aseg_gdf_attribute, netcdf_attribute in self.settings['variable_attributes'].items():
                        attribute_value = variable_attributes.get(netcdf_attribute)
                        if attribute_value is not None:
                            optional_attribute_list.append('{aseg_gdf_attribute}={attribute_value}'.format(aseg_gdf_attribute=aseg_gdf_attribute,
                                                                                                           attribute_value=attribute_value
                                                                                                           ))
                        
                    long_name = variable_attributes.get('long_name')
                    if long_name:
                        optional_attribute_list.append('NAME={long_name}'.format(long_name=long_name))
                        
                    if optional_attribute_list:
                        definition = ' , '.join(optional_attribute_list)
                    else:
                        definition = None
    
                    self.write_record2dfn_file(dfn_file,
                                               rt='',
                                               name=field_name,
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
                wkt = None
                for crs_variable_name in ['crs', 'transverse_mercator']:
                    if crs_variable_name in self.nc_dataset.variables.keys():
                        wkt = self.nc_dataset.variables[crs_variable_name].spatial_ref
                        break
                assert wkt, 'No Coordinate Reference System defined'
                
                spatial_ref = get_spatial_ref_from_wkt(wkt)
                geogcs = spatial_ref.GetAttrValue('geogcs') # e.g. 'GDA94'
                projcs = spatial_ref.GetAttrValue('projcs') # e.g. 'UTM Zone 54, Southern Hemisphere'
                ellipse_name = spatial_ref.GetAttrValue('spheroid', 0)
                major_axis = float(spatial_ref.GetAttrValue('spheroid', 1))
                prime_meridian = float(spatial_ref.GetAttrValue('primem', 1))
                inverse_flattening = float(spatial_ref.GetInvFlattening())
                #eccentricity = spatial_ref.GetAttrValue('spheroid', 2) # Non-standard definition same as inverse_flattening?
                
                if spatial_ref.IsProjected():
                    if projcs.startswith(geogcs):
                        projection_name = projcs
                    else:
                        projection_name = geogcs + ' / ' + re.sub('[\:\,\=]+', '', projcs) # e.g. 'GDA94 / UTM Zone 54, Southern Hemisphere'
                    projection_method = spatial_ref.GetAttrValue('projection').replace('_', ' ')
                    projection_parameters = [(key, float(value))
                                              for key, value in re.findall('PARAMETER\["(.+)",(\d+\.?\d*)\]', spatial_ref.ExportToPrettyWkt())
                                              ]
                else: # Unprojected CRS
                    projection_name = None
                    projection_method = None
                    projection_parameters = None

                self.defn = 0  # reset DEFN number
                
                # write 'DEFN 1 ST=RECD,RT=PROJ; RT:A4'
                self.write_record2dfn_file(dfn_file,
                                           rt='PROJ',
                                           name='RT',
                                           aseg_gdf_format='A4'
                                           )
                
                if projcs:
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
                                                   aseg_gdf_format='D14.0',
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
            logger.info('Finished writing .dfn file {}'.format(self.dfn_out_path))
        
        
        def write_dat_file():
            '''
            Helper function to output .dat file
            '''
            def line_generator():
                '''
                Generator to yield all line strings across all point variables
                '''                    
                def chunk_line_generator(start_row, end_row):
                    '''
                    Helper Generator to yield line strings for specified rows across all point variables
                    '''
                    def expand_indexing_variable(indexing_field_name, start_row, end_row):
                        '''
                        Helper function to expand indexing variables and return an array of the required size
                        '''
                        row_range = end_row - start_row
                        value_variable = self.nc_dataset.variables[indexing_field_name]
                        start_variable = self.nc_dataset.variables['index_' + indexing_field_name]
                        count_variable = self.nc_dataset.variables['index_count_' + indexing_field_name]
                        
                        expanded_array = np.zeros(shape=(row_range,), 
                                                  dtype=value_variable.dtype)
                        
                        # Assume monotonically increasing start indices to find all relevant indices
                        indices = np.where(np.logical_and((start_variable[:] >= start_row),
                                                          (start_variable[:] <= end_row)))[0]
                        
                        #logger.debug('indices: {}'.format(indices))
                        for index in indices:
                            start_index = max(start_variable[index]-start_row, 0)
                            end_index = min(start_index+count_variable[index], row_range)
                            expanded_array[start_index:end_index] = value_variable[index]
                            
                        #logger.debug('expanded_array: {}'.format(expanded_array))
                        return expanded_array
                    
                    row_range = end_row - start_row
                    column_start = 0                    
                    for field_definition in self.field_definitions:
                        variable_name = field_definition['short_name']
                        variable = self.nc_dataset.variables[variable_name]
                        
                        if len(variable.shape) == 1: # 1D variable
                            # Not an indexing variable
                            if ('point' in variable.dimensions) and (variable_name not in self.settings['index_fields']):
                                data_array = variable[start_row:end_row]
                                
                                # Include fill_values if array is masked
                                if type(data_array) == np.ma.core.MaskedArray:
                                    data_array = data_array.data
                                    
                                memory_cache_array[0:row_range, 
                                                   column_start] = variable[start_row:end_row]
                            # Indexing variable
                            elif ('point' not in variable.dimensions) and (variable_name in self.settings['index_fields']): 
                                memory_cache_array[0:row_range, 
                                                   column_start] = expand_indexing_variable(variable_name, start_row, end_row)
                            else:
                                raise BaseException('Invalid dimension for variable {}'.format(variable_name))  
                              
                        elif len(variable.shape) == 2: # 2D variable
                            data_array = variable[start_row:end_row]
                            
                            # Include fill_values if array is masked
                            if type(data_array) == np.ma.core.MaskedArray:
                                data_array = data_array.data
                                
                            memory_cache_array[0:row_range, 
                                                column_start:column_start+field_definition['columns']] = data_array
                        
                        column_start += field_definition['columns']
                         
                    for line_index in range(row_range):
                        #yield ' '.join(['{}'.format(value) for value in memory_cache_array[line_index,:]]) 
                        yield ' '.join([column_format_list[column_index].format(memory_cache_array[line_index, column_index]) 
                                        for column_index in range(total_columns)
                                        ]
                                       ) + '\n'
                        
                        
                
                # Define formats for individual columns
                column_format_list = []
                for field_definition in self.field_definitions:
                    for _column_index in range(field_definition['columns']):
                        column_format_list.append(field_definition['python_format'])
                    
                
                
                total_columns = sum([field_definition['columns']
                                     for field_definition in self.field_definitions
                                     ]
                                    )
                logger.debug('total_columns: {}'.format(total_columns))
                                
                max_chunk_size = max([field_definition['chunk_size']
                                      for field_definition in self.field_definitions
                                      ]
                                     )
                
                # Calculate maximum number of rows which can be read into memory with float64 cells
                read_chunk_size = (MAX_MEMORY_BYTES // total_columns // 8 // max_chunk_size) * max_chunk_size
                
                logger.debug('read_chunk_size: {}'.format(read_chunk_size))
                
                memory_cache_array = np.zeros(shape=(read_chunk_size, total_columns), dtype='float64')
                
                # Process all complete chunks
                point_count = 0
                for chunk_index in range(self.total_points // read_chunk_size):
                    logger.debug('Reading chunk {} for rows {}-{}'.format(chunk_index+1,
                                                                          chunk_index*read_chunk_size,
                                                                          (chunk_index+1)*read_chunk_size-1)
                                                                          )
                    for line in chunk_line_generator(chunk_index*read_chunk_size,
                                                     (chunk_index+1)*read_chunk_size):
                        point_count += 1
                        
                        if point_count == point_count // 10000 * 10000:
                            logger.info('{} points written'.format(point_count))
                            #logger.debug('line: {}'.format(line))
                    
                        yield line
                        
                        if POINT_LIMIT and (point_count >= POINT_LIMIT):
                            break
                        
                    if POINT_LIMIT and (point_count >= POINT_LIMIT):
                        break
                    
                # All complete chunks processed - process any remaining rows
                remaining_points = self.total_points % read_chunk_size
                if remaining_points and (not POINT_LIMIT or (point_count < POINT_LIMIT)):
                    logger.debug('Reading final chunk with {} points'.format(remaining_points))
                    for line in chunk_line_generator(self.total_points-remaining_points,
                                                     self.total_points):
                        point_count += 1
                        
                        if point_count == point_count // 10000 * 10000:
                            logger.info('{} points written'.format(point_count))
                            #logger.debug('line: {}'.format(line))

                        yield line
                        
                        if POINT_LIMIT and (point_count >= POINT_LIMIT):
                            break
                             
                logger.info('{} lines output'.format(point_count))
        
            # Create, write and close .dat file
            dat_file = open(self.dat_out_path, mode='w')
            for line in line_generator():
                dat_file.write(line)
            dat_file.close()
            logger.info('Finished writing .dat file {}'.format(self.dat_out_path))
                
        
        
        write_dfn_file()
        write_dat_file()
    
       
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
                            dest="dfn_out_path")
        parser.add_argument("-s", "--settings",
                            help="Path to settings file",
                            type=str,
                            dest="settings_path")
        
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='output debug information. Default is no debug info')
        
        parser.add_argument('positional_args', 
                            nargs=argparse.REMAINDER,
                            help='<nc_in_path> [<dat_out_path>]')

        return parser.parse_args()
    
    args = get_args()
    
    # Setup Logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    base_logger = logging.getLogger('geophys_utils.netcdf_converter') # Logger for base class ToNetCDFConverter
    base_logger.setLevel(level=log_level)
    logger.setLevel(level=log_level)

    assert 1 <= len(args.positional_args) <= 2, 'Invalid number of positional arguments.\n\
Usage: python {} <options> <nc_in_path> [<dat_out_path>]'.format(os.path.basename(sys.argv[0]))

    nc_in_path = args.positional_args[0]

    if len(args.positional_args) >= 2:
        dat_out_path = args.positional_args[1]
    else:
        dat_out_path = os.path.splitext(nc_in_path)[0] + '.dat'
        
    dfn_out_path = args.dfn_out_path or os.path.splitext(dat_out_path)[0] + '.dfn'
    
    netcdf2aseg_gdf_converter = NetCDF2ASEGGDFConverter(nc_in_path,
                                                        dat_out_path,
                                                        dfn_out_path
                                                        )
    
    netcdf2aseg_gdf_converter.convert2aseg_gdf()
    
       
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
        