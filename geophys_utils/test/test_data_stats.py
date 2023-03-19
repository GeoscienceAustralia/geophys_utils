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
Unit tests for geophys_utils._data_stats against a NetCDF file

Created on 15/11/2016

@author: Alex Ip
"""
import os
import unittest

from geophys_utils._data_stats import DataStats


class TestDataStats(unittest.TestCase):
    """Unit tests for geophys_utils._data_stats module."""

    NC_PATH = 'test_grid.nc'
    MAX_BYTES = 1600
    MAX_ERROR = 0.000001
    EXPECTED_RESULT = {'y_size': 178,
                       'data_type': 'float32',
                       'min': -308.20889,
                       'max': 3412.063,
                       # 'nc_path': u'C:\\Users\\u76345\\git\\geophys_utils\\geophys_utils\\test\\test_grid.nc',
                       'x_size': 79,
                       'mean': -15.958283494050137,
                       'nodata_value': -99999.0}

    def test_data_stats(self):
        print('Testing DataStats class')
        TestDataStats.EXPECTED_RESULT['nc_path'] = os.path.join(os.path.dirname(__file__), TestDataStats.NC_PATH)
        data_stats = DataStats(TestDataStats.EXPECTED_RESULT['nc_path'],
                               max_bytes=TestDataStats.MAX_BYTES)
        for key in sorted(TestDataStats.EXPECTED_RESULT.keys()):
            try:
                error = data_stats.value(key) - TestDataStats.EXPECTED_RESULT[key]
                assert error <= TestDataStats.MAX_ERROR, 'Incorrect numeric value for {}\nExpected {}, got {}.'.format(
                    key,
                    TestDataStats.EXPECTED_RESULT[key],
                    data_stats.value(key))
            except TypeError:
                assert data_stats.value(key) == TestDataStats.EXPECTED_RESULT[key], 'Incorrect value for %s' % key


# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestDataStats]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())


if __name__ == '__main__':
    main()
