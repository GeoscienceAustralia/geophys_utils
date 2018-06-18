'''
Created on 13 Jun. 2018

@author: u76345
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
from math import log10, ceil

from geophys_utils.netcdf_converter import ToNetCDFConverter, NetCDFVariable
from geophys_utils import get_spatial_ref_from_wkt

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

# Maximum size of in-memory cache array (in bytes)
MAX_MEMORY_BYTES = 1073741824 # 1GB

# Default effective chunk size for un-chunked variables.
# Set to 0 for no chunking (i.e. complete read)
DEFAULT_READ_CHUNK_SIZE = 1
    
class NetCDF2ASEGGDFConverter(object):
    '''
    NetCDF2ASEGGDFConverter class definition to convert netCDF file to ASEG-GDF format
    '''
    # Approximate number of significant figures for each datatype
    SIG_FIGS = OrderedDict([('int8', 2),
                            ('int16', 4),
                            ('int32', 10),
                            ('int64', 19),
                            ('float32', 6),
                            ('float64', 20)
                            ]
                           )

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
        
        self.settings_path = settings_path or os.path.splitext(__file__)[0] + '_settings.yml'
        
        try:
            self.settings = yaml.safe_load(open(self.settings_path))
        except:
            self.settings = {}
            
        logger.debug('self.settings: {}'.format(pformat(self.settings)))

        self.nc_dataset = netCDF4.Dataset(self.netcdf_in_path, 'r')
        assert 'point' in self.nc_dataset.dimensions.keys(), 'No point dimension defined in netCDF dataset'
        
        # Build field definitions
        self.field_definitions = OrderedDict()
        for variable_name, variable in self.nc_dataset.variables.items():
            
            # Skip scalar variables or those which don't have a point dimension
            if not len(variable.dimensions) or ('point' not in variable.dimensions):
                continue
            
            chunk_size = variable.chunking()[0]
            # If variable is not chunked and DEFAULT_READ_CHUNK_SIZE is defined
            if DEFAULT_READ_CHUNK_SIZE and (chunk_size == variable.shape[0]):
                chunk_size = min(DEFAULT_READ_CHUNK_SIZE, variable.shape[0])
                
            if len(variable.shape) == 1: # 1D variable
                columns = 1
            elif len(variable.shape) == 2: # 2D variable
                columns = variable.shape[1]
                
            dtype = str(variable.dtype)
            
            sig_figs = NetCDF2ASEGGDFConverter.SIG_FIGS[dtype] # Look up approximate significant figures
            integer_digits = ceil(log10(np.nanmax(variable[:]) + 1.0))
            fractional_digits = sig_figs - integer_digits
            
            if dtype.startswith('int'):
                fmt = 'I{}'.format(integer_digits)
            elif dtype.startswith('float'):
                fmt = 'F{}.{}'.format(integer_digits, fractional_digits)
            #===================================================================
            # elif dtype.startswith('float'):
            #     fmt = 'E{}.{}'.format(integer_digits, fractional_digits)
            #===================================================================
            
            # Pre-pend column count to start of fmt
            if columns > 1:
                fmt = '{}{}'.format(columns, fmt)
            
            #TODO: Add extra field definition stuff like ASEG-GDF format specifier
            field_definition = {'dtype': str(variable.dtype),
                                'chunk_size': chunk_size,
                                'columns': columns,
                                'fmt': fmt,
                                'integer_digits': integer_digits,
                                'fractional_digits': fractional_digits
                                }
            
            self.field_definitions[variable_name] = field_definition

        logger.debug('self.field_definitions: {}'.format(pformat(self.field_definitions)))
    
    def convert2aseg_gdf(self): 
        '''
        Function to convert netCDF file to ASEG-GDF
        '''
        def write_dfn_file():
            '''
            Helper function to output .dfn file
            '''
            # Create, write and close .dat file
            dfn_file = open(self.dfn_out_path, 'w')
            dfn_file.write('DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76\n') # TODO: Check this first line 
            
            defn = 0
            for field_name, field_definition in self.field_definitions.items():
                defn += 1
                
                variable = self.nc_dataset.variables[field_name]
                #print(field_name, variable.dtype)
                
                line = 'DEFN {defn} ST=RECD,RT=; {short_name} : {fmt}'.format(defn=defn,
                                                                              short_name=field_name,
                                                                              fmt=field_definition['fmt']
                                                                              )
                
                variable_attributes = variable.__dict__
                
                optional_attribute_list = []
                
                units = variable_attributes.get('units')
                if units:
                    optional_attribute_list.append('UNITS = {units}'.format(units=units))

                long_name = variable_attributes.get('long_name')
                if long_name:
                    optional_attribute_list.append(long_name)
                    
                if optional_attribute_list:
                    line += ' : ' + ' , '.join(optional_attribute_list) 

                dfn_file.write(line + '\n')
                
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
                    row_range = end_row - start_row
                    column_start = 0                    
                    for field_name, field_definition in self.field_definitions.items():
                        variable = self.nc_dataset.variables[field_name]
                        
                        if len(variable.shape) == 1: # 1D variable
                            memory_cache_array[0:row_range, 
                                               column_start] = variable[start_row:end_row]
                        if len(variable.shape) == 2: # 2D variable
                            memory_cache_array[0:row_range, 
                                                column_start:column_start+field_definition['columns']] = variable[start_row:end_row]
                        column_start += field_definition['columns']
                         
                    for line_index in range(row_range):
                        #TODO: Fix format - temporary implementation with comma-separated strings
                        yield ','.join(['{}'.format(value) for value in memory_cache_array[line_index,:]]) 
                
                point_count = self.nc_dataset.dimensions['point'].size
                
                total_columns = sum([field_definition['columns']
                                     for field_definition in self.field_definitions.values()
                                     ]
                                    )
                logger.debug('total_columns: {}'.format(total_columns))
                                
                max_chunk_size = max([field_definition['chunk_size']
                                      for field_definition in self.field_definitions.values()
                                      ]
                                     )
                
                # Calculate maximum number of rows which can be read into memory with float64 cells
                read_chunk_size = (MAX_MEMORY_BYTES // total_columns // 8 // max_chunk_size) * max_chunk_size
                
                logger.debug('read_chunk_size: {}'.format(read_chunk_size))
                
                memory_cache_array = np.zeros(shape=(read_chunk_size, total_columns), dtype='float64')
                
                # Process all complete chunks
                for chunk_index in range(point_count // read_chunk_size):
                    logger.debug('Reading chunk {} for rows {}-{}'.format(chunk_index+1,
                                                                          chunk_index*read_chunk_size,
                                                                          (chunk_index+1)*read_chunk_size-1)
                                                                          )
                    for line in chunk_line_generator(chunk_index*read_chunk_size,
                                                     (chunk_index+1)*read_chunk_size):
                        yield line
                        
                # All complete chunks processed - process any remaining rows
                remaining_points = point_count % read_chunk_size
                if remaining_points:
                    logger.debug('Reading final chunk with {} points'.format(remaining_points))
                    for line in chunk_line_generator(point_count-remaining_points,
                                                     point_count):
                        yield line
                             
                logger.info('{} lines output'.format(self.nc_dataset.dimensions['point'].size))
        
            # Create, write and close .dat file
            dat_file = open(self.dat_out_path, 'w')
            for line in line_generator():
                dat_file.write(line + '\n')
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
        