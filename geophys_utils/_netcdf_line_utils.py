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
import os
import numpy as np
from geophys_utils._netcdf_point_utils import NetCDFPointUtils
from scipy.spatial.distance import pdist
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
                 enable_disk_cache=None,
                 enable_memory_cache=True,
                 cache_dir=None,
                 debug=False):
        '''
        NetCDFLineUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a line dataset
        @parameter enable_disk_cache: Boolean parameter indicating whether local cache file should be used, or None for default 
        @parameter enable_memory_cache: Boolean parameter indicating whether values should be cached in memory or not.
        @parameter debug: Boolean parameter indicating whether debug output should be turned on or not
        '''     
        # Start of init function - Call inherited constructor first
        NetCDFPointUtils.__init__(self, netcdf_dataset=netcdf_dataset, 
                                  enable_disk_cache=enable_disk_cache, 
                                  enable_memory_cache=enable_memory_cache,
                                  cache_dir=cache_dir,
                                  debug=debug)

        self._line = None
        self._line_index = None
            
        
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
                line_index = int(np.where(self.line == line_number)[0])
            except TypeError:
                logger.warning('Invalid line number %d' % line_number)
                continue # Should this be break?
            
            line_mask = (self.line_index == line_index)
            #logger.debug('Line {} has a total of {} points'.format(line_number, np.count_nonzero(line_mask))) 
            
            yield line_number, line_mask
    
    
    def get_lines(self, line_numbers=None, 
                  variables=None, 
                  bounds=None, 
                  bounds_wkt=None,
                  subsampling_distance=None
                  ):
        '''
        Generator to return coordinates and specified variable values for specified lines
        @param line_numbers: list of integer line number or single integer line number
        @param variables: list of variable name strings or single variable name string. None returns all variables
        @param bounds: Spatial bounds for point selection
        @param bounds_wkt: WKT for bounds Coordinate Reference System 
        @param subsampling_distance: Minimum subsampling_distance expressed in native coordinate units (e.g. degrees)
        
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

        # Return all variables if specified variable is None
        variables = self.point_variables if variables is None else variables
        
        # Allow single variable to be given as a string
        single_var = (type(variables) == str)
        if single_var:
            variables = [variables]
        
        bounds = bounds or self.bounds
        
        spatial_subset_mask = self.get_spatial_mask(self.get_reprojected_bounds(bounds, bounds_wkt, self.wkt))
        
        logger.debug('subsampling_distance: {}'.format(subsampling_distance))
        for line_number, line_mask in self.get_line_masks(line_numbers=line_numbers):
        
            point_indices = np.where(np.logical_and(line_mask, spatial_subset_mask))[0]
            #logger.debug('Line {} has {} points in bounding box'.format(line_number, len(point_indices))) 
            line_point_count = len(point_indices)
            if line_point_count:
                # Use subset of indices if stride is set
                if subsampling_distance:
                    line_length = pdist([self.xycoords[point_indices[0]], self.xycoords[point_indices[-1]]])[0]
                    logger.debug('line_length: {}'.format(line_length))
                    stride = int(line_point_count/max(1, line_length/subsampling_distance))
                    logger.debug('stride: {}'.format(stride))
                    # Create array of subset indices, including the index of the last point if not already in subsample indices
                    subset_indices = np.unique(np.concatenate((np.arange(0, line_point_count, stride),
                                                              np.array([line_point_count-1])),
                                                              axis=None)
                                                              )
                    logger.debug('Subset of line {} has {} points'.format(line_number, len(subset_indices)))
                    point_indices = point_indices[subset_indices]
                    
                line_dict = {'coordinates': self.xycoords[point_indices]}
                # Add <variable_name>: <variable_array> for each specified variable
                for variable_name in variables:
                    line_dict[variable_name] = self.netcdf_dataset.variables[variable_name][point_indices]
        
                yield line_number, line_dict
 
    
    def get_line_values(self):
        '''
        Function to retrieve array of line number values from self.netcdf_dataset
        '''
        line_variable = self.netcdf_dataset.variables.get('line')
        assert line_variable, 'Variable "line" does not exist in netCDF file'
        if line_variable.shape: # Multiple lines
            #line_values = self.fetch_array(line_variable)
            line_values = line_variable[:] # Should be small enough to retrieve in one hit
            
        else: # Scalar - only one line
            line_values = line_variable[:].reshape((1,)) # Change scalar to single-element 1D array
            
        return line_values
    
    def get_line_index_values(self):
        '''
        Function to retrieve array of line_index indices from self.netcdf_dataset
        '''
        if self.line.shape: # Multiple lines
            line_index_variable = self.netcdf_dataset.variables.get('line_index')
            if line_index_variable: # Lookup format lines - Current format
                logger.debug('Line data is in lookup format (current)')
                #line_indices = self.fetch_array(line_index_variable)
                line_indices = line_index_variable[:]
            else: # Indexing format lines - OLD FORMAT
                raise BaseException('Line data is in indexing format (unsupported)')
            
        else: # Scalar
            # Synthesize line_indices array with all zeroes for single value
            line_indices = np.zeros(shape=(self.point_count,), 
                                    dtype='int8'
                                    ) 
            
        return line_indices
      
    
    @property
    def line(self):
        '''
        Property getter function to return array of all line numbers
        Always cache this in memory - should only be small
        '''
        if self._line is None:
            self._line = self.get_line_values()
            
        return self._line

    @property
    def line_index(self):
        '''
        Property getter function to return line_indices for all points
        '''
        if self.enable_memory_cache and self._line_index is not None:
            #logger.debug('Returning memory cached line_index')
            return self._line_index
            
        if self.enable_disk_cache:
            line_index_path = self.cache_basename + '_lineindex.dat'
            
            if os.path.isfile(line_index_path):
                # Cached line_index file exists - read it
                logger.debug('Reading line_index cache file {}'.format(line_index_path))
                with open(line_index_path, 'rb') as line_index_file:
                    line_index = np.fromfile(line_index_file, dtype=np.int16)
            else:
                line_index = self.get_line_index_values()
                logger.debug('Saving line_index cache file {}'.format(line_index_path))
                os.makedirs(self.cache_dir, exist_ok=True)
                with open(line_index_path, 'wb') as line_index_file:
                    line_index.astype(np.int16).tofile(line_index_file) # Write to cache file
            
        else: # No caching - read line_index from source file
            line_index = self.get_line_index_values()

        if self.enable_memory_cache:
            self._line_index = line_index
            
        return line_index
    
