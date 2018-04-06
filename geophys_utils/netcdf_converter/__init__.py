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
NetCDFConverter abstract base class for converting data to netCDF

Created on 28Mar.2018

@author: Alex Ip
'''
import abc
import logging
import netCDF4
import numpy as np
import osgeo
from collections import OrderedDict

from geophys_utils import get_spatial_ref_from_wkt

class NetCDFVariable(object):
    '''
    Class to manage netCDF variable contents
    '''
    DEFAULT_VARIABLE_PARAMETERS = {'complevel': 4, 
                           'zlib': True, 
                           'fletcher32': True,
                           'shuffle': True,
                           'endian': 'little',
                           }
    
    def __init__(self, short_name, data, dimensions, fill_value, attributes, dtype=None):
        '''
        Constructor for class NetCDFVariable to manage netCDF variable contents
        @param dtype: Optional datatype to override data
        '''
        self.short_name = short_name # String used for variable name
        self.data = data # Numpy array
        self.dimensions = dimensions # List of <dimension_name> strings for array, or None or empty list for scalar
        self.fill_value = fill_value # Fill value for no-data array values
        self.attributes = attributes # dict of variable attribute <key>:<value> pairs 
        self.dtype = dtype or (data.dtype if (type(data) == np.ndarray) else 'float64')
        #TODO: Implement something better for determining the type of scalars
                
        
    def create_var_in_dataset(self, nc_output_dataset, chunk_size=None):
        '''
        Function to create netCDF variable in specified dataset
        '''
        variable_parameters = NetCDFVariable.DEFAULT_VARIABLE_PARAMETERS
        
        # If not a scalar, check array shape against specified dimensions in netCDF dataset
        if self.dimensions:
            assert set(self.dimensions) <= set(list(nc_output_dataset.dimensions.keys())), 'Invalid dimension(s) specified'

            assert self.data.shape == tuple([nc_output_dataset.dimensions[dimension_name].size
                                             for dimension_name in self.dimensions
                                             ]), 'Invalid array shape for specified dimension(s)'
                                             
            if chunk_size:
                variable_parameters['chunksizes'] = [min(nc_output_dataset.dimensions[dimension_name].size,
                                                         chunk_size
                                                         )
                                                    for dimension_name in self.dimensions
                                                    ]
                
            if self.fill_value is not None:
                variable_parameters['fill_value'] = self.fill_value
                
                                             
        output_variable = nc_output_dataset.createVariable(self.short_name,
                                                           self.dtype,
                                                           self.dimensions,
                                                           **variable_parameters
                                                           )
        # Set data values
        output_variable[:] = self.data
        
        # Set variable attributes
        output_variable.setncatts(self.attributes)


class NetCDFConverter(object):
    '''
    NetCDFConverter abstract base class for converting data to netCDF
    '''
    # Define single default chunk size for all dimensions
    DEFAULT_CHUNK_SIZE = 1024
    
    @abc.abstractmethod
    def __init__(self, nc_out_path, netcdf_format='NETCDF4_CLASSIC', chunk_size=None):
        '''
        Abstract base constructor for abstract base class NetCDFConverter
        Needs to initialise object with everything that is required for the other abstract base methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        @param nc_out_path: Path to output netCDF file on filesystem
        @param netcdf_format: Format for netCDF file. Defaults to 'NETCDF4_CLASSIC'
        @param chunk_size: single default chunk size for all dimensions
        '''
        self.nc_out_path = nc_out_path
        
        self.chunk_size = chunk_size or NetCDFConverter.DEFAULT_CHUNK_SIZE
        
        # Create netCDF output file
        self.nc_output_dataset = netCDF4.Dataset(nc_out_path, 
                                                 mode="w", 
                                                 clobber=True, 
                                                 format=netcdf_format)
        
    
    def __del__(self):
        '''
        Concrete destructor for abstract base class NetCDFConverter
        '''
        self.nc_output_dataset.close()
        
    @abc.abstractmethod
    def get_global_attributes(self):
        '''
        Abstract base method to return dict of global attribute <key>:<value> pairs       
        '''
        return {}
    
    @abc.abstractmethod
    def get_dimensions(self):
        '''
        Abstract base method to return OrderedDict of <dimension_name>:<dimension_size> pairs       
        '''
        dimensions = OrderedDict()
        
        #=======================================================================
        # # Example lat/lon dimensions
        # dimensions['lon'] = 509
        # dimensions['lat'] = 639  
        #=======================================================================
              
        return dimensions
    
    @abc.abstractmethod
    def variable_generator(self):
        '''
        Abstract base generator to yield NetCDFVariable objects       
        '''        
        #=======================================================================
        # # Example of latitude dimension variable creation
        # yield self.build_dim_index_variable(dimension_name='lat', 
        #                                     min_value=-22.9247209891964, 
        #                                     max_value=-20.5641209891964, 
        #                                     long_name='latitude', 
        #                                     units='degrees north', 
        #                                     standard_name='latitude',
        #                                     descending=True # Invert Y axis
        #                                     )
        # 
        # # Example of longitude dimension variable creation
        # yield self.build_dim_index_variable(dimension_name='lon', 
        #                                     min_value=121.122089060582, 
        #                                     max_value=123.001689060582, 
        #                                     long_name='longitude', 
        #                                     units='degrees east', 
        #                                     standard_name='longitude',
        #                                     descending=False
        #                                     )
        #=======================================================================

#===============================================================================
#         # Example of crs variable creation for GDA94
#         yield self.build_crs_variable('''\
# GEOGCS["GDA94",
#     DATUM["Geocentric_Datum_of_Australia_1994",
#         SPHEROID["GRS 1980",6378137,298.257222101,
#             AUTHORITY["EPSG","7019"]],
#         TOWGS84[0,0,0,0,0,0,0],
#         AUTHORITY["EPSG","6283"]],
#     PRIMEM["Greenwich",0,
#         AUTHORITY["EPSG","8901"]],
#     UNIT["degree",0.0174532925199433,
#         AUTHORITY["EPSG","9122"]],
#     AUTHORITY["EPSG","4283"]]
# '''
#             )
#===============================================================================
        return
    
    
    def preprocess_netcdf(self):
        '''
        Function to perform any pre-processing on netCDF file before dimensions and variables
        have been created. May be overridden in subclass.
        '''
        return
    
    def postprocess_netcdf(self):
        '''
        Function to perform any post-processing on netCDF file after dimensions and variables
        have been created. May be overridden in subclass.
        '''
        return
    
    def build_dim_index_variable(self, 
                                 dimension_name, 
                                 min_value, 
                                 max_value, 
                                 long_name, 
                                 units, 
                                 standard_name=None,
                                 descending=False
                                 ):
        '''
        Concrete method to build dimension index variable
        N.B: Need to have dimension defined prior to calling this
        
        @param dimension_name: Name of dimension 
        @param min_value: Minimum value in dimension variable array
        @param max_value: Maximum value in dimension variable array
        @param long_name: long_name attribute string
        @param units: units attribute string
        @param standard_name: Optional standard_name attribute string
        @param descending: Boolean flag determining whether array values are descending. Default = False
        '''
        #=======================================================================
        # # Example dimension variables
        # double lat(lat) ;
        #         lat:units = "degrees_north" ;
        #         lat:long_name = "latitude" ;
        #         lat:standard_name = "latitude" ;
        # double lon(lon) ;
        #         lon:units = "degrees_east" ;
        #         lon:long_name = "longitude" ;
        #         lon:standard_name = "longitude" ;
        #=======================================================================
        assert dimension_name in self.nc_output_dataset.dimensions.keys(), 'Invalid dimension name'
        
        dimension_variable_attributes = {'long_name': long_name,
                                         'units': units
                                         }
        
        if standard_name:
            dimension_variable_attributes['standard_name'] = standard_name
        
        increment = (max_value - min_value) / self.nc_output_dataset.dimensions[dimension_name].size
        
        if descending:
            data_array = np.arange(max_value, min_value, -increment, dtype='float64')
        else:
            data_array = np.arange(min_value, max_value, increment, dtype='float64')
        
        return NetCDFVariable(short_name=dimension_name, 
                              data=data_array, 
                              dimensions=[dimension_name], 
                              fill_value=None, 
                              attributes=dimension_variable_attributes)    
    
    def build_crs_variable(self, crs, grid_dimensions=None):
        '''
        Concrete method to build "crs" or "transverse_mercator" NetCDFVariable instance from well known text
        N.B: Need to create dimensions and dimension variables first if grid_dimensions is specified
        @param crs: Either osgeo.osr.SpatialReference or WKT string defining Coordinate Reference System
        @grid_dimensions: list of two dimension names for spatial grid, or None for ungridded data
        @grid_resolutions: list of two floats defining resolutions for spatial grid, or None for ungridded data
        '''
        def get_geotransform(grid_dimensions):
            '''
            Helper function to return geotransform for gridded data.
            Assumes yx array ordering for grid
            N.B: This will fail if dimensions and dimension variables don't exist
            '''
            assert len(grid_dimensions) == 2, 'grid_dimensions must be of length 2'
            assert set(grid_dimensions) <= set(self.nc_output_dataset.dimensions.keys()), 'Invalid grid_dimensions'
            assert set(grid_dimensions) <= set(self.nc_output_dataset.variables.keys()), 'Dimension index variables not created'
            
            cell_sizes = [(self.nc_output_dataset.variables[grid_dimensions[dim_index]][-1] -
                           self.nc_output_dataset.variables[grid_dimensions[dim_index]][0]) /
                           self.nc_output_dataset.dimensions[grid_dimensions[dim_index]].size
                           for dim_index in range(2)
                          ]
            
            return [self.nc_output_dataset.variables[grid_dimensions[1]][0] - cell_sizes[1] / 2, # x_min
                    cell_sizes[1], # x_size
                    0.0,
                    self.nc_output_dataset.variables[grid_dimensions[0]][0] - cell_sizes[1] / 2, # y_min
                    0.0,
                    cell_sizes[0], # y_size
                    ]
            
            
        # Determine wkt and spatial_ref from crs as required
        if type(crs) == osgeo.osr.SpatialReference:
            spatial_ref = crs
        elif type(crs) == str:
            spatial_ref = get_spatial_ref_from_wkt(crs)
            
        wkt = spatial_ref.ExportToWkt() # Export WKT from spatial_ref for consistency even if supplied
        crs_attributes = {'spatial_ref': wkt}
        
        crs_attributes['inverse_flattening'] = spatial_ref.GetInvFlattening()
        crs_attributes['semi_major_axis'] = spatial_ref.GetSemiMajor()
        crs_attributes['longitude_of_prime_meridian'] = spatial_ref.GetAttrValue('PRIMEM', 1)

        if spatial_ref.GetUTMZone(): # CRS is UTM
            crs_variable_name = 'transverse_mercator'
            crs_attributes['grid_mapping_name'] = 'transverse_mercator'

            crs_attributes['latitude_of_projection_origin'] = spatial_ref.GetProjParm('latitude_of_origin')
            crs_attributes['scale_factor_at_central_meridian'] = spatial_ref.GetProjParm('scale_factor')
            crs_attributes['longitude_of_central_meridian'] = spatial_ref.GetProjParm('central_meridian')
            crs_attributes['false_northing'] = spatial_ref.GetProjParm('false_northing')
            crs_attributes['false_easting'] = spatial_ref.GetProjParm('false_easting')
                   
        else:
            #===================================================================
            # # Example crs attributes created by GDAL:
            #
            # crs:inverse_flattening = 298.257222101 ;
            # crs:spatial_ref = "GEOGCS[\"GEOCENTRIC DATUM of AUSTRALIA\",DATUM[\"GDA94\",SPHEROID[\"GRS80\",6378137,298.257222101]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]" ;
            # crs:semi_major_axis = 6378137. ;
            # crs:GeoTransform = "121.1202390605822 0.0037 0 -20.56227098919639 0 -0.0037 " ;
            # crs:grid_mapping_name = "latitude_longitude" ;
            # crs:longitude_of_prime_meridian = 0. ;
            #===================================================================
            crs_variable_name = 'crs'
            crs_attributes['grid_mapping_name'] = 'latitude_longitude'
            
        # Set GeoTransform only for regular gridded data
        if grid_dimensions:
            crs_attributes['GeoTransform'] = ' '.join([str(value) for value in get_geotransform(grid_dimensions)]) 

        return NetCDFVariable(short_name=crs_variable_name, 
                              data=0, 
                              dimensions=[], # Scalar
                              fill_value=None, 
                              attributes=crs_attributes,
                              dtype='int8' # Byte datatype
                              )    

    def convert2netcdf(self):
        '''
        Concrete method to output netCDF
        '''
        self.preprocess_netcdf()
        
        # Create dimensions in netCDF output file
        for dimension_name, dimension_size in iter(self.get_dimensions().items()):
            self.nc_output_dataset.createDimension(dimname=dimension_name, size=dimension_size)
            
        # Create variables in netCDF output file
        for nc_variable in self.variable_generator():
            nc_variable.create_var_in_dataset(self.nc_output_dataset, chunk_size=self.chunk_size)
            
        # Set global attributes in netCDF output file
        for attribute_name, attribute_value in iter(self.get_global_attributes().items()):
            setattr(self.nc_output_dataset, attribute_name, attribute_value)
            
        self.postprocess_netcdf()
            
            