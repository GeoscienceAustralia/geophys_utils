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
"""
Unit tests for geophys_utils._netcdf_line_utils against a NetCDF line data file

Created on 15/11/2016

@author: Alex Ip
"""
import os
import re
import unittest

import netCDF4
import numpy as np
from shapely.geometry.polygon import Polygon

from geophys_utils._netcdf_line_utils import NetCDFLineUtils

netcdf_line_utils = None

NC_PATH = 'https://dapds00.nci.org.au/thredds/dodsC/iv65/Geoscience_Australia_Geophysics_Reference_Data_Collection/airborne_geophysics/SA/line/P1255/P1255-line-magnetic-Marree-AWAGS_MAG_2010.nc'
NC_TITLE = 'Marree, SA, 2012 (P1255), magnetic line data, AWAGS levelled'

TEST_BOUNDS = (137, -29, 138, -28)

MAX_BYTES = 1600
MAX_ERROR = 0.000001

TEST_GET_LINE_MASK_RESULTS = ((190520, 5032),
                              (190500, 4994)
                              )

TEST_GET_LINE_RESULTS = ((190520, 7, 10064),
                         (190500, 7, 9988)
                         )


class TestNetCDFLineUtilsConstructor(unittest.TestCase):
    """Unit tests for TestNetCDFLineUtils Constructor.
    N.B: This should be run first"""

    def test_netcdf_line_utils_constructor(self):
        print('Testing NetCDFLineUtils constructor')
        global netcdf_line_utils

        if re.match('^http.*', NC_PATH):
            nc_path = NC_PATH
        else:
            nc_path = os.path.join(os.path.dirname(__file__), NC_PATH)
        print(nc_path)
        nc_dataset = netCDF4.Dataset(nc_path)
        netcdf_line_utils = NetCDFLineUtils(nc_dataset)

        # print(netcdf_line_utils.__dict__)
        assert nc_dataset.title == NC_TITLE, 'Invalid dataset title: "{}" != "{}"'.format(nc_dataset.title, NC_TITLE)


class TestNetCDFLineUtilsFunctions1(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""

    def test_get_line_masks(self):
        print('Testing get_line_masks function')
        count = 0
        for line_number, line_mask in netcdf_line_utils.get_line_masks():
            # print('Line {} has {} points'.format(line_number, np.count_nonzero(line_mask)))
            assert (line_number, np.count_nonzero(line_mask)) == TEST_GET_LINE_MASK_RESULTS[
                count], "Invalid get_line_masks result"
            count += 1
            if count >= 2:
                break

    def test_concave_hull(self):
        print('Testing concave hull')
        assert isinstance(netcdf_line_utils.get_concave_hull(), Polygon)


class TestNetCDFLineUtilsFunctions2(unittest.TestCase):
    """Unit tests for geophys_utils._netcdf_line_utils functions"""

    def test_get_lines(self):
        print('Testing get_lines function')
        count = 0
        for line_number, line_dict in netcdf_line_utils.get_lines():
            # ===================================================================
            # print('Line {} has {} variables with {} points'.format(line_number,
            #                                                    len(line_dict)-1, 
            #                                                    np.count_nonzero(line_dict['coordinates'])
            #                                                    )
            #       )
            # ===================================================================
            assert (line_number, len(line_dict) - 1, np.count_nonzero(line_dict['coordinates'])) == \
                   TEST_GET_LINE_RESULTS[count], \
                "Invalid get_lines result: Expected {}, got {}".format(TEST_GET_LINE_RESULTS[count], (
                line_number, len(line_dict) - 1, np.count_nonzero(line_dict['coordinates'])))

            count += 1
            if count >= 2:
                break


# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestNetCDFLineUtilsConstructor,
                    TestNetCDFLineUtilsFunctions1,
                    TestNetCDFLineUtilsFunctions2
                    ]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())


if __name__ == '__main__':
    main()
