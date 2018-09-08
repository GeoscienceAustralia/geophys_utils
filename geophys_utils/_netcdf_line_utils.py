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
Created on 16/11/2016

@author: Alex Ip
'''
import numpy as np
from geophys_utils._netcdf_point_utils import NetCDFPointUtils
import logging

# Setup logging handlers if required
logger = logging.getLogger(__name__) # Get __main__ logger
logger.setLevel(logging.INFO) # Initial logging level for this module
    
class NetCDFLineUtils(NetCDFPointUtils):
    '''
    NetCDFLineUtils class to do various fiddly things with NetCDF geophysics line data files.
    '''

    def __init__(self, 
                 netcdf_dataset, 
                 debug=False):
        '''
        NetCDFLineUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a line dataset
        '''     
        # Start of init function - Call inherited constructor first
        #TODO: Make local caching optional - currently forced to True
        NetCDFPointUtils.__init__(self, netcdf_dataset=netcdf_dataset, 
                                  enable_cache=True, 
                                  debug=debug)

        line_variable = self.netcdf_dataset.variables.get('line')
        assert line_variable, 'Variable "line" does not exist in netCDF file'
        if line_variable.shape: # Multiple lines
            line_values = self.fetch_array(line_variable)
            self.line_count = line_variable.shape[0]
            line_index_variable = self.netcdf_dataset.variables.get('line_index')
            if line_index_variable: # Lookup format lines - Current format
                logger.debug('Line data is in lookup format (current)')
                line_indices = self.fetch_array(line_index_variable)
                line_index_dtype = line_index_variable.dtype
                line_index_var_options = line_index_variable.filters() or {}
            else: # Indexing format lines - OLD FORMAT
                #TODO: Get rid of this case when we stop using old format data
                index_line_variable = self.netcdf_dataset.variables.get('index_line')
                index_count_variable = self.netcdf_dataset.variables.get('index_count')
                assert (index_line_variable.dimensions[0] == 'line'
                        and index_count_variable.dimensions[0] == 'line'), 'Invalid line variable dimensioning'
                logger.warning('Line data is in indexing format (deprecated)')
                
                # Synthesise line index array from start & count values
                line_indices = np.zeros(shape=(self.point_count,), 
                                        dtype=('int8' if len(line_values) <= 127
                                            else 'int16' if len(line_values) <= 32768
                                                else 'int32')
                                        )
                for line_index in range(len(line_values)):
                    line_indices[index_line_variable[line_index]:index_line_variable[line_index]+index_count_variable[line_index]] = line_index         
                
                line_index_dtype = index_line_variable.dtype
                line_index_var_options = index_line_variable.filters() or {}
        else: # Scalar
            line_values = line_variable[:].reshape((1,)) # Change scalar to 1D array
            self.line_count = 1
            # Synthesize line_indices array with all zeroes for single value
            line_indices = np.zeros(shape=(self.point_count,), 
                                    dtype=('int8' if len(line_values) <= 127
                                            else 'int16' if len(line_values) <= 32768
                                                else 'int32')
                                    ) 
            line_index_dtype = 'int8'
            line_index_var_options = {}
        
        self._nc_cache_dataset.createDimension('line', self.line_count if not self.unlimited_points else None)
        
        line_var_options = line_variable.filters() or {}
        line_var_options['zlib'] = False
        # if hasattr(self.netcdf_dataset.variables['line'], '_FillValue'):
        #     line_var_options['fill_value'] = self.netcdf_dataset.variables['line']._FillValue
            
        line_cache_variable = self._nc_cache_dataset.createVariable('line', 
                                      line_variable.dtype, 
                                      ('line'),
                                      **line_var_options
                                      )
        self.line = line_cache_variable
        self.line[:] = line_values
        
        line_index_var_options['zlib'] = True
        # if hasattr(self.netcdf_dataset.variables['line'], '_FillValue'):
        #     line_var_options['fill_value'] = self.netcdf_dataset.variables['line']._FillValue
            
        line_index_cache_variable = self._nc_cache_dataset.createVariable('line_index', 
                                      line_index_dtype, 
                                      ('point'),
                                      **line_index_var_options
                                      )
        self.line_index = line_index_cache_variable
        self.line_index[:] = line_indices
        
        self.kdtree = None
        
        print(line_values)
        print(line_indices)
        
    def get_line_masks(self, line_numbers=None):
        '''
        Generator to return boolean masks of dimension 'point' for specified lines
        @param line_numbers: list of integer line number or single integer line number
        
        @return line_number: line number for single line
        @return line_mask: boolean mask for single line

        '''
        # Yield masks for all lines if not specified
        if line_numbers is None:
            line_numbers = self.line[:]

        # Convert single line number to single element list
        try:
            _line_numbers_iterator = iter(line_numbers)
        except TypeError:
            line_numbers = [line_numbers]

        for line_number in line_numbers:
            try:
                line_index = int(np.where(self.line[:] == line_number)[0])
            except TypeError:
                logger.warning('Invalid line number %d' % line_number)
                continue # Should this be break?
            
            line_mask = self.line_index[:] == line_index
            logger.debug('line_mask: {}'.format(line_mask))  
            
            yield line_number, line_mask
    
    
    def get_lines(self, line_numbers=None, variables=None, bounds=None, bounds_wkt=None):
        '''
        Generator to return coordinates and specified variable values for specified lines
        @param line_numbers: list of integer line number or single integer line number
        @param variables: list of variable name strings or single variable name string
        @param bounds: Spatial bounds for point selection
        @param bounds_wkt: WKT for bounds Coordinate Reference System 
        
        @return line_number: line number for single line
        @return: dict containing coords and values for required variables keyed by variable name
        '''
        # Return all lines if not specified
        if line_numbers is None:
            line_numbers = self.line[:]

        # Convert single line number to single element list
        try:
            _line_numbers_iterator = iter(line_numbers)
        except TypeError:
            line_numbers = [line_numbers]

        # Allow single variable to be given as a string
        variables = variables or self.point_variables
        single_var = (type(variables) == str)
        if single_var:
            variables = [variables]
        
        bounds = bounds or self.bounds
        
        spatial_subset_mask = self.get_spatial_mask(self.get_reprojected_bounds(bounds, bounds_wkt, self.wkt))
        
        for line_number in line_numbers:
            _line_number, line_mask = next(self.get_line_masks(line_numbers=line_number)) # Only one mask per line
        
            point_indices = np.where(np.logical_and(line_mask, spatial_subset_mask))[0]
            if len(point_indices):
                line_dict = {'coordinates': self.xycoords[point_indices]}
                for variable_name in variables:
                    line_dict[variable_name] = self.netcdf_dataset.variables[variable_name][point_indices]
        
                yield line_number, line_dict

