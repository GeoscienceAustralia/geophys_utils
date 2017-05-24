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
from distutils.core import setup
import os

version = '0.0.0'

setup(name='geophys_utils',
      version=version,
      packages=[
          'geophys_utils',
          'geophys_utils.test'
      ],
      package_data={'geophys_utils': ['csw_utils_settings.yml']
                    },
      scripts=(['bin/csw_find',
                'bin/rechunk'] 
               if (os.name == 'posix')
               else (['bin\\csw_find.bat',
                      'bin\\rechunk.bat'] 
                     if (os.name == 'nt')
                     else [])),
      requires=[
            'distutils',
            'functools',
            'itertools',
            'netCDF4',
            'numpy',
            'osgeo',
            'owslib',
            'scipy',
            'shapely',
            'tempfile',
            'unittest',
            'yaml'
            ],
      url='https://github.com/alex-ip/geophys_utils',
      author='Alex Ip - Geoscience Australia',
      maintainer='Alex Ip - Geoscience Australia',
      maintainer_email='alex.ip@ga.gov.au',
      description='Geophysics data access utilities',
      long_description='Geophysics data access utilities',
      license='Creative Commons Attribution 4.0 International'
      )
