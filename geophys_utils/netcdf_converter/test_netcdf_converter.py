#!/usr/bin/env python

# ===============================================================================
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
# ===============================================================================
'''
TestNetCDFConverter concrete class for converting data to netCDF

Created on 28Mar.2018

@author: Alex Ip
'''
from collections import OrderedDict

import numpy as np

from geophys_utils.netcdf_converter import ToNetCDFConverter, NetCDFVariable


class TestNetCDFConverter(ToNetCDFConverter):
    '''
    TestNetCDFConverter concrete class for converting CSV data to netCDF
    '''

    def __init__(self, nc_out_path, netcdf_format='NETCDF4_CLASSIC'):
        '''
        Concrete constructor for subclass TestNetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure the base class constructor is called from the subclass constructor
        '''
        ToNetCDFConverter.__init__(self, nc_out_path, netcdf_format)

    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs       
        '''
        return {'title': 'test dataset'}

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs       
        '''
        dimensions = OrderedDict()

        # Example lat/lon dimensions
        dimensions['lon'] = 509
        dimensions['lat'] = 639

        return dimensions

    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects       
        '''
        # Example of latitude dimension variable creation
        yield self.build_dim_index_variable(dimension_name='lat',
                                            min_value=-22.9247209891964,
                                            max_value=-20.5641209891964,
                                            long_name='latitude',
                                            units='degrees north',
                                            standard_name='latitude',
                                            descending=True  # Invert Y axis
                                            )

        # Example of longitude dimension variable creation
        yield self.build_dim_index_variable(dimension_name='lon',
                                            min_value=121.122089060582,
                                            max_value=123.001689060582,
                                            long_name='longitude',
                                            units='degrees east',
                                            standard_name='longitude',
                                            descending=False
                                            )

        # Example of crs variable creation for GDA94
        yield self.build_crs_variable('''\
GEOGCS["GDA94",
    DATUM["Geocentric_Datum_of_Australia_1994",
        SPHEROID["GRS 1980",6378137,298.257222101,
            AUTHORITY["EPSG","7019"]],
        TOWGS84[0,0,0,0,0,0,0],
        AUTHORITY["EPSG","6283"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.0174532925199433,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4283"]]
'''
                                      )

        yield NetCDFVariable(short_name='test_data',
                             data=np.random.random((self.nc_output_dataset.dimensions['lat'].size,
                                                    self.nc_output_dataset.dimensions['lon'].size)),
                             dimensions=['lat', 'lon'],
                             fill_value=0.0,
                             attributes={'units': 'random crap',
                                         'long_name': 'random numbers between 0 and 1'
                                         },
                             dtype='float32'
                             )

        return


def main():
    nc_out_path = 'C:\\Temp\\test.nc'
    tnc = TestNetCDFConverter(nc_out_path)
    tnc.convert2netcdf()
    print('Finished writing netCDF file {}'.format(nc_out_path))


if __name__ == '__main__':
    main()
