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

    def __init__(self, netcdf_dataset, debug=False):
        '''
        NetCDFLineUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a line dataset
        '''     
        # Start of init function - Call inherited constructor first
        NetCDFPointUtils.__init__(self, netcdf_dataset, debug)

        line_dimension = self.netcdf_dataset.dimensions['line']
        self.line_count = len(line_dimension)
        
        self._nc_cache_dataset.createDimension('line', self.line_count if not self.unlimited_points else None)
        
        var_options = self.netcdf_dataset.variables['line'].filters() or {}
        var_options['zlib'] = False
        if hasattr(self.netcdf_dataset.variables['line'], '_FillValue'):
            var_options['fill_value'] = self.netcdf_dataset.variables['line']._FillValue
            
        self._nc_cache_dataset.createVariable('line', 
                                      self.netcdf_dataset.variables['line'].dtype, 
                                      ('line'),
                                      **var_options
                                      )
        self.line = self._nc_cache_dataset.variables['line']
        self.line[:] = self.fetch_array(self.netcdf_dataset.variables['line'])
        
        self._nc_cache_dataset.createDimension('start_end', 2)
        
        # Cater for two possible variable naming conventions
        #TODO: Remove this after the variable naming scheme has been finalised
        try:
            index_line_variable = self.netcdf_dataset.variables['index_line']
            index_count_variable = self.netcdf_dataset.variables['index_count']
        except:
            index_line_variable = self.netcdf_dataset.variables['_index_line']
            index_count_variable = self.netcdf_dataset.variables['_index_count']
        
        var_options = index_line_variable.filters() or {}
        
        var_options['zlib'] = False
        if hasattr(index_line_variable, '_FillValue'):
            var_options['fill_value'] = index_line_variable._FillValue
            
        self._nc_cache_dataset.createVariable('line_start_end', 
                                      index_line_variable.dtype, 
                                      ('line', 'start_end'),
                                      **var_options
                                      )
        self.line_start_end = self._nc_cache_dataset.variables['line_start_end']
        self.line_start_end[:,0] = self.fetch_array(index_line_variable)
        self.line_start_end[:,1] = self.fetch_array(index_count_variable)
        self.line_start_end[:,1] += self.line_start_end[:,0]
        
        self.kdtree = None
        
    def get_line_masks(self, line_numbers=None):
        '''
        Generator to return boolean masks of dimension 'point' for specified lines
        @param line_numbers: list of integer line number or single integer line number
        
        @return line_number: line number for single line
        @return line_mask: boolean mask for single line

        '''
        line_number_array = self.line[...]
        line_start_end_array = self.line_start_end[...]
        
        # Yield masks for all lines if not specified
        if line_numbers is None:
            line_numbers = line_number_array

        # Convert single line number to single element list
        try:
            _line_numbers_iterator = iter(line_numbers)
        except TypeError:
            line_numbers = [line_numbers]

        for line_number in line_numbers:
            try:
                line_index = int(np.where(line_number_array == line_number)[0])
            except TypeError:
                logger.warning('Invalid line number %d' % line_number)
                continue # Should this be break?
            
            line_mask = np.zeros((self.point_count,), bool)
            line_mask[line_start_end_array[line_index,0]:line_start_end_array[line_index,1]] = True
            
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
            line_numbers = self.line[...]

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
