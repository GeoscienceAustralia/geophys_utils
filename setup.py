#!/usr/bin/env python

from distutils.core import setup
import os

version = '0.0.0'

setup(name='geophys_utils',
      version=version,
      packages=[
          'geophys_utils',
          'geophys_utils.test'
      ],
      package_data={
      },
      scripts=(['bin/csw_find'] if (os.name == 'posix')
               else (['bin\\csw_find.bat'] if (os.name == 'nt')
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
            ],
      url='https://github.com/alex-ip/geophys_utils',
      author='Alex Ip - Geoscience Australia',
      maintainer='Alex Ip - Geoscience Australia',
      maintainer_email='alex.ip@ga.gov.au',
      description='Geophysics data access utilities',
      long_description='Geophysics data access utilities',
      license='Creative Commons Attribution 4.0 International'
      )
