#!/usr/bin/env python

from distutils.core import setup

version = '0.0.0'

setup(name='geophys_utils',
      version=version,
      packages=[
          'geophys_utils',
          'geophys_utils.test'
      ],
      package_data={
      },
      scripts=[
          ],
      requires=[
            'distutils',
            'functools',
            'itertools',
            'logging',
            'math',
            'netCDF4',
            'numpy',
            'os',
            'osgeo',
            'owslib',
            're',
            'scipy',
            'shapely',
            'sys',
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
