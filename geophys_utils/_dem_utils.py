#!/usr/bin/env python

#===============================================================================
# Copyright (c)  2014 Geoscience Australia
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither Geoscience Australia nor the names of its contributors may be
#       used to endorse or promote products derived from this software
#       without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#===============================================================================

'''
Created on 05/10/2012

@author: Alex Ip
'''
import os
import sys
import logging
import numpy
import gc
from osgeo import osr
import numexpr
import netCDF4
import gc
from scipy.ndimage import sobel
from geophys_utils._netcdf_grid_utils import NetCDFGridUtils 
from geophys_utils._array_pieces import array_pieces 

from _vincenty import vinc_dist
from _blrb import interpolate_grid
    
# Set top level standard output 
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(message)s')
console_handler.setFormatter(console_formatter)

logger = logging.getLogger(__name__)
if not logger.level:
    logger.setLevel(logging.DEBUG) # Default logging level for all modules
    logger.addHandler(console_handler)
                
RADIANS_PER_DEGREE = 0.01745329251994329576923690768489
    
class earth(object):

    # Mean radius
    RADIUS = 6371009.0  # (metres)

    # WGS-84
    #RADIUS = 6378135.0  # equatorial (metres)
    #RADIUS = 6356752.0  # polar (metres)

    # Length of Earth ellipsoid semi-major axis (metres)
    SEMI_MAJOR_AXIS = 6378137.0

    # WGS-84
    A = 6378137.0           # equatorial radius (metres)
    B = 6356752.3142        # polar radius (metres)
    F = (A - B) / A         # flattening
    ECC2 = 1.0 - B**2/A**2  # squared eccentricity

    MEAN_RADIUS = (A*2 + B) / 3

    # Earth ellipsoid eccentricity (dimensionless)
    #ECCENTRICITY = 0.00669438
    #ECC2 = math.pow(ECCENTRICITY, 2)

    # Earth rotational angular velocity (radians/sec)
    OMEGA = 0.000072722052

class DEMUtils(NetCDFGridUtils):

    def getFileSizekB(self, path):
        """Gets the size of a file (megabytes).
    
        Arguments:
            path: file path
     
        Returns:
            File size (MB)
    
        Raises:
            OSError [Errno=2] if file does not exist
        """     
        return os.path.getsize(path) / 1024

    def getFileSizeMB(self, path):
        """Gets the size of a file (megabytes).
    
        Arguments:
            path: file path
     
        Returns:
            File size (MB)
    
        Raises:
            OSError [Errno=2] if file does not exist
        """     
        return self.getFileSizekB(path) / 1024
    
    def get_pixel_size(self, (x, y)):
        """
        Returns X & Y sizes in metres of specified pixel as a tuple.
        N.B: Pixel ordinates are zero-based from top left
        """
        logger.debug('(x, y) = (%f, %f)', x, y)
        
        native_spatial_reference = osr.SpatialReference()
        native_spatial_reference.ImportFromWkt(self.crs)
        
        latlong_spatial_reference = native_spatial_reference.CloneGeogCS()
        coord_transform_to_latlong = osr.CoordinateTransformation(native_spatial_reference, latlong_spatial_reference)

        # Determine pixel centre and edges in georeferenced coordinates
        xw = self.GeoTransform[0] + x * self.GeoTransform[1]
        yn = self.GeoTransform[3] + y * self.GeoTransform[5] 
        xc = self.GeoTransform[0] + (x + 0.5) * self.GeoTransform[1]
        yc = self.GeoTransform[3] + (y + 0.5) * self.GeoTransform[5] 
        xe = self.GeoTransform[0] + (x + 1.0) * self.GeoTransform[1]
        ys = self.GeoTransform[3] + (y + 1.0) * self.GeoTransform[5] 
        
        logger.debug('xw = %f, yn = %f, xc = %f, yc = %f, xe = %f, ys = %f', xw, yn, xc, yc, xe, ys)
        
        # Convert georeferenced coordinates to lat/lon for Vincenty
        lon1, lat1, _z = coord_transform_to_latlong.TransformPoint(xw, yc, 0)
        lon2, lat2, _z = coord_transform_to_latlong.TransformPoint(xe, yc, 0)
        logger.debug('For X size: (lon1, lat1) = (%f, %f), (lon2, lat2) = (%f, %f)', lon1, lat1, lon2, lat2)
        x_size, _az_to, _az_from = vinc_dist(earth.F, earth.A, 
                                             lat1 * RADIANS_PER_DEGREE, lon1 * RADIANS_PER_DEGREE, 
                                             lat2 * RADIANS_PER_DEGREE, lon2 * RADIANS_PER_DEGREE)
        
        lon1, lat1, _z = coord_transform_to_latlong.TransformPoint(xc, yn, 0)
        lon2, lat2, _z = coord_transform_to_latlong.TransformPoint(xc, ys, 0)
        logger.debug('For Y size: (lon1, lat1) = (%f, %f), (lon2, lat2) = (%f, %f)', lon1, lat1, lon2, lat2)
        y_size, _az_to, _az_from = vinc_dist(earth.F, earth.A, 
                                             lat1 * RADIANS_PER_DEGREE, lon1 * RADIANS_PER_DEGREE, 
                                             lat2 * RADIANS_PER_DEGREE, lon2 * RADIANS_PER_DEGREE)
        
        logger.debug('(x_size, y_size) = (%f, %f)', x_size, y_size)
        return (x_size, y_size)

    def get_pixel_size_grid(self, source_array, offsets, dimension_index):
        """ Returns grid with interpolated X or Y pixel sizes for given datasets"""
        
        def get_pixel_size(x, y):
            return self.get_pixel_size((offsets[0] + x, offsets[1] + y))[dimension_index]
        
        pixel_size_grid = numpy.zeros(source_array.shape).astype(source_array.dtype)
        interpolate_grid(depth=1, shape=pixel_size_grid.shape, eval_func=get_pixel_size, grid=pixel_size_grid)
        
        return pixel_size_grid


    def __init__(self, dem_dataset):
        """Constructor
        Arguments:
            source_dem_nc: NetCDF dataset containing DEM data
        """    
        # Start of init function - Call inherited constructor first
        NetCDFGridUtils.__init__(self, dem_dataset)        

    def create_dzdxy_arrays(self, elevation_array, offsets):
        '''
        Function to return two arrays containing dzdx and dzdy values
        '''
        def pixels_in_m(): 
            ''' 
            Function returning True if pixels are in metres
            '''
            result = True
            for dimension_name in self.data_variable.dimensions:
                try:
                    if self.netcdf_dataset.variables[dimension_name].units == 'm':
                        continue
                    else:
                        result = False
                        break
                except:
                    result = False
                    break
                
            return result
                    
        native_pixel_x_size = float(abs(self.GeoTransform[1]))
        native_pixel_y_size = float(abs(self.GeoTransform[5]))

        dzdx_array = sobel(elevation_array, axis=1)/(8. * native_pixel_x_size)
        dzdy_array = sobel(elevation_array, axis=0)/(8. * native_pixel_y_size)

        if pixels_in_m():  
            print 'Pixels are a uniform size of %f x %f metres.' % (native_pixel_x_size, native_pixel_y_size)       
            # Pixel sizes are in metres - use scalars
            pixel_x_metres = native_pixel_x_size
            pixel_y_metres = native_pixel_y_size
        else:  
            print 'Pixels are of varying sizes. Computing and applying pixel size arrays.'       
            # Compute variable pixel size  
            m_array = self.get_pixel_size_grid(elevation_array, offsets)
            pixel_x_metres = m_array[:,:,0]
            pixel_y_metres = m_array[:,:,1]
         
            dzdx_array = numexpr.evaluate("dzdx_array * native_pixel_x_size / pixel_x_metres")
            dzdy_array = numexpr.evaluate("dzdy_array * native_pixel_y_size / pixel_y_metres")
        
        return dzdx_array, dzdy_array
    
    def create_slope_array(self, dzdx_array, dzdy_array):
        hypotenuse_array = numpy.hypot(dzdx_array, dzdy_array)
        slope_array = numexpr.evaluate("arctan(hypotenuse_array) / RADIANS_PER_DEGREE")
        #Blank out no-data cells
        slope_array[numpy.isnan(slope_array)] = self.data_variable._FillValue
        
        return slope_array
    
        
    def create_aspect_array(self, dzdx_array, dzdy_array):
        # Convert angles from conventional radians to compass heading 0-360
        aspect_array = numexpr.evaluate("(450 - arctan2(dzdy_array, -dzdx_array) / RADIANS_PER_DEGREE) % 360")
        #Blank out no-data cells
        aspect_array[numpy.isnan(aspect_array)] = self.data_variable._FillValue
        
        return aspect_array

        
    def create_slope_and_aspect(self, slope_path=None, aspect_path=None, overlap=2):
        '''
        Create slope & aspect datasets from elevation
        '''
        # Copy dataset structure but not data
        slope_path = slope_path or os.path.splitext(self.nc_path)[0] + '_slope.nc'
        self.copy(slope_path, empty_var_list=[self.data_variable.name])
        slope_nc_dataset = netCDF4.Dataset(slope_path, 'r+')
        slope_nc_dataset.renameVariable(self.data_variable.name, 'slope')
        slope_variable = slope_nc_dataset.variables['slope']
        slope_variable.long_name = 'slope expressed in degrees from horizontal (0=horizontal, 90=vertical)'
        slope_variable.units = 'degrees'

        aspect_path = aspect_path or os.path.splitext(self.nc_path)[0] + '_aspect.nc'
        self.copy(aspect_path, empty_var_list=[self.data_variable.name])
        aspect_nc_dataset = netCDF4.Dataset(aspect_path, 'r+')
        aspect_nc_dataset.renameVariable(self.data_variable.name, 'aspect')
        aspect_variable = aspect_nc_dataset.variables['aspect']
        aspect_variable.long_name = 'aspect expressed compass bearing of normal to plane (0=North, 90=East, etc.)'
        aspect_variable.units = 'degrees'
              
        # Process dataset in small pieces
        for piece_array, offsets in array_pieces(self.data_variable, 
                                                 max_bytes=self.max_bytes if self.opendap else self.max_bytes/2, # Need to allow for multiple arrays in memory 
                                                 overlap=overlap):
            print 'Processing array of shape %s at %s' % (piece_array.shape, offsets)
            
            if type(piece_array) == numpy.ma.masked_array:
                piece_array = piece_array.data # Convert from masked array to plain array

            piece_array[(piece_array == self.data_variable._FillValue)] = numpy.NaN
            
            # Calculate raw source & destination slices including overlaps
            source_slices = [slice(0, 
                                   piece_array.shape[dim_index])
                             for dim_index in range(2)
                             ]
                                  
            dest_slices = [slice(offsets[dim_index], 
                                 offsets[dim_index] + piece_array.shape[dim_index])
                           for dim_index in range(2)
                           ]
                                  
            # Trim overlaps off source & destination slices
            source_slices = [slice(0 if dest_slices[dim_index].start < overlap else source_slices[dim_index].start+overlap, 
                                   piece_array.shape[dim_index] if (self.data_variable.shape[dim_index] - dest_slices[dim_index].stop) < overlap else source_slices[dim_index].stop-overlap)
                             for dim_index in range(2)
                             ]
                                  
            dest_slices = [slice(0 if dest_slices[dim_index].start < overlap else dest_slices[dim_index].start+overlap, 
                                 self.data_variable.shape[dim_index] if (self.data_variable.shape[dim_index] - dest_slices[dim_index].stop) < overlap else dest_slices[dim_index].stop-overlap)
                           for dim_index in range(2)
                           ]
            
            print 'Computing dzdx and dzdy arrays'
            dzdx_array, dzdy_array = self.create_dzdxy_arrays(piece_array, offsets)
            
            print 'Computing slope array'
            result_array = self.create_slope_array(dzdx_array, dzdy_array)

            print 'Writing slope array of shape %s at %s' % (tuple([dest_slices[dim_index].stop - dest_slices[dim_index].start
                                                     for dim_index in range(2)
                                                     ]),
                                              tuple([dest_slices[dim_index].start
                                                     for dim_index in range(2)
                                                     ])
                                              )
            slope_variable[dest_slices] = result_array[source_slices]  
            slope_nc_dataset.sync()
                 
            print 'Computing aspect array'
            result_array = self.create_aspect_array(dzdx_array, dzdy_array)
                                  
            print 'Writing aspect array of shape %s at %s' % (tuple([dest_slices[dim_index].stop - dest_slices[dim_index].start
                                                     for dim_index in range(2)
                                                     ]),
                                              tuple([dest_slices[dim_index].start
                                                     for dim_index in range(2)
                                                     ])
                                              )
            aspect_variable[dest_slices] = result_array[source_slices]      
            aspect_nc_dataset.sync()   
            
        slope_nc_dataset.close() 
        print 'Finished writing slope dataset %s' % slope_path
        
        aspect_nc_dataset.close()     
        print 'Finished writing aspect dataset %s' % aspect_path
                
if __name__ == '__main__':
    # Define command line arguments
    dem_path = sys.argv[1] 
    
    try:
        slope_path = sys.argv[2] 
    except:
        slope_path = None
        
    try:
        aspect_path = sys.argv[3]
    except:
        aspect_path = None
    
    dem_utils = DEMUtils(dem_path)
    
    dem_utils.create_slope_and_aspect(slope_path, aspect_path)
