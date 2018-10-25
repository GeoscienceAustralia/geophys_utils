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
import netCDF4

# Setup logging handlers if required
logger = logging.getLogger(__name__) # Get __main__ logger
logger.setLevel(logging.INFO) # Initial logging level for this module
    
class NetCDFLineUtils(NetCDFPointUtils):
    '''
    NetCDFLineUtils class to do various fiddly things with NetCDF geophysics line data files.
    '''

    def __init__(self, 
                 netcdf_dataset,
                 memcached_connection=None,
                 enable_disk_cache=None,
                 enable_memory_cache=True,
                 cache_path=None,
                 debug=False):
        '''
        NetCDFLineUtils Constructor
        @parameter netcdf_dataset: netCDF4.Dataset object containing a line dataset
        @parameter enable_disk_cache: Boolean parameter indicating whether local cache file should be used, or None for default 
        @parameter enable_memory_cache: Boolean parameter indicating whether values should be cached in memory or not.
        @parameter debug: Boolean parameter indicating whether debug output should be turned on or not
        '''     
        # Start of init function - Call inherited constructor first
        super().__init__(netcdf_dataset=netcdf_dataset, 
                         memcached_connection=memcached_connection,
                         enable_disk_cache=enable_disk_cache, 
                         enable_memory_cache=enable_memory_cache,
                         cache_path=cache_path,
                         debug=debug)

        logger.debug('Running NetCDFLineUtils constructor')

        # Initialise private property variables to None until set by property getter methods
        self._line = None
        self._line_index = None
            
        
    def get_line_masks(self, line_numbers=None, subset_mask=None, get_contiguous_lines=False):
        '''
        Generator to return boolean masks of dimension 'point' for specified lines
        @param line_numbers: list of integer line number or single integer line number, or None for all lines
        @param subset_mask: optional Boolean mask for subset (e.g. spatial mask)
        @param get_contiguous_lines: Boolean flag indicating whether masked gaps in lines should be included
        
        @return line_number: line number for single line
        @return line_mask: Boolean mask for single line

        '''       
        # Yield masks for all lines in subset if no line numbers specified
        if line_numbers is None:
            line_number_subset = self.line # All line numbers
        else:
            # Convert single line number to single element list
            try:
                _line_numbers_iterator = iter(line_numbers)
            except TypeError:
                line_numbers = [line_numbers]

            line_number_subset = np.array(line_numbers)
            
        if subset_mask is not None:
            line_number_subset = line_number_subset[np.isin(line_number_subset, self.line[np.unique(self.line_index[subset_mask])])] # Exclude lines not in subset
        else:    
            line_number_subset = line_number_subset[np.isin(line_number_subset, self.line)] # Exclude bad line numbers 
    
        line_mask = np.zeros(shape=self.line_index.shape, dtype=np.bool) # Keep re-using same in-memory array

        for line_number in line_number_subset:
            line_mask[:] = False
            line_index = int(np.where(self.line == line_number)[0])
            
            if subset_mask is not None:
                line_mask[subset_mask] = (self.line_index[subset_mask] == line_index)
                
                if get_contiguous_lines:
                    # Include all points in line from first to last in subset
                    point_indices = np.where(line_mask)[0]
                    line_mask[min(point_indices):max(point_indices)+1] = True
            else:
                line_mask[(self.line_index == line_index)] = True
                
            #logger.debug('Line {} has a total of {} points'.format(line_number, np.count_nonzero(line_mask))) 
            
            if np.any(line_mask): # This is probably redundant
                yield line_number, line_mask
    
    
    def get_lines(self, line_numbers=None, 
                  variables=None, 
                  bounds=None, 
                  #bounds_wkt=None,
                  subsampling_distance=None,
                  get_contiguous_lines=False
                  ):
        '''
        Generator to return coordinates and specified variable values for specified lines
        @param line_numbers: list of integer line number or single integer line number
        @param variables: list of variable name strings or single variable name string. None returns all variables
        @param bounds: Spatial bounds for point selection
        @param bounds_wkt: WKT for bounds Coordinate Reference System 
        @param subsampling_distance: Minimum subsampling_distance expressed in native coordinate units (e.g. degrees)
        @param get_contiguous_lines: Boolean flag indicating whether masked gaps in lines should be included
        
        @return line_number: line number for single line
        @return: dict containing coords and values for required variables keyed by variable name
        '''
        # Return all variables if specified variable is None
        variables = self.point_variables if variables is None else variables
        
        # Allow single variable to be given as a string
        single_var = (type(variables) == str)
        if single_var:
            variables = [variables]
        
        bounds = bounds or self.bounds
        
        #spatial_subset_mask = self.get_spatial_mask(self.get_reprojected_bounds(bounds, bounds_wkt, self.wkt))
        spatial_subset_mask = self.get_spatial_mask(bounds)
        
        logger.debug('subsampling_distance: {}'.format(subsampling_distance))
        
        for line_number, line_mask in self.get_line_masks(line_numbers=line_numbers, 
                                                          subset_mask=spatial_subset_mask,
                                                          get_contiguous_lines=get_contiguous_lines
                                                          ):
        
            point_indices = np.where(line_mask)[0]
            #logger.debug('Line {} has {} points in bounding box'.format(line_number, len(point_indices))) 
            line_point_count = len(point_indices)
            if line_point_count: # This test should be redundant
                # Use subset of indices if stride is set
                if subsampling_distance:
                    line_length = pdist([self.xycoords[point_indices[0]], self.xycoords[point_indices[-1]]])[0]
                    logger.debug('line_length: {}'.format(line_length))
                    stride = max(1, int(line_point_count/max(1, line_length/subsampling_distance)))
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
        if len(self.netcdf_dataset.variables['line']): # Multiple lines
            line_index_variable = self.netcdf_dataset.variables.get('line_index')
            if line_index_variable: # Lookup format lines - Current format
                line_indices = self.fetch_array(line_index_variable)
                #line_indices = line_index_variable[:]
            else: # Indexing format lines - OLD FORMAT
                raise BaseException('Line data is in indexing format (unsupported)')
            
        else: # Scalar
            # Synthesize line_indices array with all zeroes for single value
            line_indices = np.zeros(shape=(self.point_count,), 
                                    dtype='int8'
                                    ) 
            
        return line_indices
      
    def get_cached_line_arrays(self):
        '''
        Helper function to cache both line & line_index
        '''
        line = None
        line_index = None

        if self.enable_disk_cache:
            if os.path.isfile(self.cache_path):
                # Cached coordinate file exists - read it
                cache_dataset = netCDF4.Dataset(self.cache_path, 'r')
                
                #assert cache_dataset.source == self.nc_path, 'Source mismatch: cache {} vs. dataset {}'.format(cache_dataset.source, self.nc_path)
                
                if 'line' in cache_dataset.variables.keys():
                    line = cache_dataset.variables['line'][:]
                    logger.debug('Read {} lines from cache file {}'.format(line.shape[0], self.cache_path))
                else:
                    logger.debug('Unable to read line variable from netCDF cache file {}'.format(self.cache_path))                

                if 'line_index' in cache_dataset.variables.keys():
                    line_index = cache_dataset.variables['line_index'][:]
                    logger.debug('Read {} line_indices from cache file {}'.format(line_index.shape[0], self.cache_path))
                else:
                    logger.debug('Unable to read line variable from netCDF cache file {}'.format(self.cache_path))  
                                  
                cache_dataset.close()
            else:
                logger.debug('NetCDF cache file {} does not exist'.format(self.cache_path))


            if line is None or line_index is None:
                if line is None:
                    line = self.get_line_values()
                if line_index is None:
                    line_index = self.get_line_index_values()
                
                os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
                if os.path.isfile(self.cache_path):
                    cache_dataset = netCDF4.Dataset(self.cache_path, 'r+')
                else:
                    cache_dataset = netCDF4.Dataset(self.cache_path, 'w')
                
                if not hasattr(cache_dataset, 'source'):
                    cache_dataset.source = self.nc_path
                    
                #assert cache_dataset.source == self.nc_path, 'Source mismatch: cache {} vs. dataset {}'.format(cache_dataset.source, self.nc_path)
                
                if 'point' not in cache_dataset.dimensions.keys():
                    cache_dataset.createDimension(dimname='point', size=line_index.shape[0]) 
                    
                if 'line' not in cache_dataset.dimensions.keys():
                    cache_dataset.createDimension(dimname='line', size=line.shape[0]) 
                    
                if 'line' not in cache_dataset.variables.keys():
                    cache_dataset.createVariable('line',
                                                 line.dtype,
                                                 dimensions=['line'],
                                                 **self.CACHE_VARIABLE_PARAMETERS
                                                 )
                cache_dataset.variables['line'][:] = line # Write lines to cache file
                
                if 'line_index' not in cache_dataset.variables.keys():
                    cache_dataset.createVariable('line_index',
                                                 line_index.dtype,
                                                 dimensions=['point'],
                                                 **self.CACHE_VARIABLE_PARAMETERS
                                                 )
                cache_dataset.variables['line_index'][:] = line_index # Write lines to cache file
                
                cache_dataset.close() 
                logger.debug('Saved {} lines for {} points to cache file {}'.format(line.shape[0], line_index.shape[0], self.cache_path))

        return line, line_index
        
    @property
    def line(self):
        '''
        Property getter function to return array of all line numbers
        Always cache this in memory - should only be small
        The order of priority for retrieval is memory, memcached, disk cache then dataset.
        '''
        line = None
        line_index = None
        if self.enable_memory_cache and self._line is not None:
            #logger.debug('Returning memory cached line')
            return self._line

        elif self.memcached_connection is not None:
            line_cache_key = self.cache_basename + '_line'
            # line = self.memcached_connection.get(line_cache_key)
            # if line is not None:
            #     logger.debug('memcached key found at {}'.format(line_cache_key))
            # else:
            #     line = self.get_line_index_values()
            #     logger.debug('Memcached key not found. Adding value with key {}'.format(line_cache_key))
            #     self.memcached_connection.add(line_cache_key, line)
            line = self.memcached_connection.get(line_cache_key)
            if line is not None:    
                logger.debug('memcached key found at {}'.format(line_cache_key))
            else:
                line = self.get_line_values()
                logger.debug('memcached key not found. Adding entry with key {}'.format(line_cache_key))
                self.memcached_connection.add(line_cache_key, line)
            
        elif self.enable_disk_cache:
            line, line_index = self.get_cached_line_arrays()           
        else: # No caching - read line from source file
            line = self.get_line_values()
            line_index = None

        if self.enable_memory_cache:
            self._line = line
            if line_index is not None:
                self._line_index = line_index
                
        #logger.debug('line: {}'.format(line))
        return line

    @property
    def line_index(self):
        '''
        Property getter function to return line_indices for all points
        The order of priority for retrieval is memory, memcached, disk cache then dataset.
        '''
        line = None
        line_index = None
        if self.enable_memory_cache and self._line_index is not None:
            #logger.debug('Returning memory cached line_index')
            return self._line_index

        elif self.memcached_connection is not None:
            line_index_cache_key = self.cache_basename + '_line_index'

            # line_index = self.memcached_connection.get(line_index_cache_key)
            # if line_index is not None:
            #     logger.debug('memcached key found at {}'.format(line_index_cache_key))
            # else:
            #     line_index = self.get_line_index_values()
            #     logger.debug('Memcached key not found. Adding value with key {}'.format(line_index_cache_key))
            #     self.memcached_connection.add(line_index_cache_key, line_index)

            line_index = self.memcached_connection.get(line_index_cache_key)
            if line_index is not None:
                logger.debug('memcached key found at {}'.format(line_index_cache_key))
            else:
                line_index = self.get_line_index_values()
                logger.debug('memcached key not found. Adding entry with key {}'.format(line_index_cache_key))
                self.memcached_connection.add(line_index_cache_key, line_index)

        elif self.enable_disk_cache:
            line, line_index = self.get_cached_line_arrays()           
        else:  # No caching - read line_index from source file
            line = None
            line_index = self.get_line_index_values()

        if self.enable_memory_cache:
            if line is not None:
                self._line = line
            self._line_index = line_index
            
        #logger.debug('line_index: {}'.format(line_index))
        return line_index
    
