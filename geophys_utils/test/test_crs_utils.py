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
"""
Unit tests for geophys_utils._crs_utils against a NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import numpy as np
import re
from osgeo.osr import CoordinateTransformation
from geophys_utils._crs_utils import get_coordinate_transformation, get_utm_wkt, transform_coords

COORD_ERROR_MARGIN = 1.0e-08

class TestCRSUtils(unittest.TestCase):
    """Unit tests for geophys_utils._crs_utils module."""
    
    EPSG4326_WKT = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'
    EPSG4326_EPSG = 'EPSG:4326'

    EPSG3577_WKT = 'PROJCS["GDA94 / Australian Albers",GEOGCS["GDA94",DATUM["Geocentric_Datum_of_Australia_1994",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6283"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4283"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["latitude_of_center",0],PARAMETER["longitude_of_center",132],PARAMETER["standard_parallel_1",-18],PARAMETER["standard_parallel_2",-36],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3577"]]'
    EPSG3577_EPSG = 'EPSG:3577'

    #UTM_WKT = 'PROJCS["UTM Zone 55,Southern Hemisphere",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["Meter",1]]'
    UTM_WKT = 'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'

    GDA2020_WKT = 'GEOGCS["GDA2020",DATUM["Geocentric_Datum_of_Australia_2020",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","1168"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","7844"]]'
    GDA2020_EPSG = 'EPSG:7844'

    EPSG4326_COORDS = [149.160, -35.306]
    #EPSG4326_COORD_ARRAY = np.array([[149.160, -35.306], [150, -35]])
    EPSG4326_COORD_ARRAY = [[149.160, -35.306], [150, -35]]
    UTM_COORDS = [696382.56321712, 6090881.85849329] #(696382.5632178178, 6090881.858493158)
    UTM_COORD_ARRAY = [(696382.5632171195, 6090881.858493287), (773798.0963396085, 6122843.308355326)]

    def test_get_coordinate_transformation(self):
        print('Testing get_coordinate_transformation function')
        #TODO: Make this one test with different test cases
        coordinate_transformation = get_coordinate_transformation(TestCRSUtils.EPSG4326_WKT,
                                                                  TestCRSUtils.EPSG3577_WKT)
        assert coordinate_transformation is not None
        assert type(coordinate_transformation) == CoordinateTransformation

        coordinate_transformation = get_coordinate_transformation(TestCRSUtils.EPSG4326_WKT,
                                                                  TestCRSUtils.GDA2020_WKT)
        assert coordinate_transformation is not None
        assert type(coordinate_transformation) == CoordinateTransformation

        coordinate_transformation = get_coordinate_transformation(TestCRSUtils.GDA2020_WKT,
                                                                  TestCRSUtils.EPSG4326_WKT
                                                                  )
        assert coordinate_transformation is not None
        assert type(coordinate_transformation) == CoordinateTransformation

        coordinate_transformation = get_coordinate_transformation(TestCRSUtils.EPSG4326_EPSG,
                                                                  TestCRSUtils.EPSG3577_EPSG)
        assert coordinate_transformation is not None
        assert type(coordinate_transformation) == CoordinateTransformation

        coordinate_transformation = get_coordinate_transformation(TestCRSUtils.EPSG4326_WKT, 
                                                                  TestCRSUtils.EPSG4326_WKT)
        assert coordinate_transformation is None, 'Null transformation should return None'

    def test_get_utm_wkt(self):
        print('Testing get_utm_wkt function')
        utm_wkt = get_utm_wkt(TestCRSUtils.EPSG4326_COORDS, 
                              TestCRSUtils.EPSG4326_EPSG)
        utm_wkt = re.sub(',\s+', ',', re.sub('\s+', ' ', utm_wkt))
        expected_wkt = re.sub(',\s+', ',', re.sub('\s+', ' ', TestCRSUtils.UTM_WKT))
        assert utm_wkt == expected_wkt, 'Incorrect UTM CRS: {} instead of {}'.format(utm_wkt, expected_wkt)

    def test_transform_coords(self):
        print('Testing transform_coords function with single coordinate{}'.format(TestCRSUtils.EPSG4326_COORDS))
        utm_coords = transform_coords(TestCRSUtils.EPSG4326_COORDS, TestCRSUtils.EPSG4326_WKT, TestCRSUtils.UTM_WKT)
        assert (abs(utm_coords - np.array(TestCRSUtils.UTM_COORDS)) < COORD_ERROR_MARGIN).all(), \
            'Incorrect UTM coordinates: {} instead of {}'.format(utm_coords, np.array(TestCRSUtils.UTM_COORDS))
        
        print('Testing transform_coords function with multi coordinate {}'.format(TestCRSUtils.EPSG4326_COORD_ARRAY))
        utm_coord_array = transform_coords(TestCRSUtils.EPSG4326_COORD_ARRAY, TestCRSUtils.EPSG4326_WKT, TestCRSUtils.UTM_WKT)
        assert (abs(utm_coord_array - TestCRSUtils.UTM_COORD_ARRAY) < COORD_ERROR_MARGIN).all(), \
            'Incorrect UTM coordinates: {} instead of {}'.format(utm_coord_array, TestCRSUtils.UTM_COORD_ARRAY)
# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestCRSUtils]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())

if __name__ == '__main__':
    main()
