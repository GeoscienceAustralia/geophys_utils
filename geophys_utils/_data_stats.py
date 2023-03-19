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
Created on 15Aug.,2016

@author: Alex
'''
import os
import sys

import netCDF4
import numpy as np

from geophys_utils._array_pieces import array_pieces


class DataStats(object):
    '''
    DataStats class definition. Obtains statistics for gridded data
    '''
    key_list = ['nc_path', 'data_type', 'nodata_value', 'x_size', 'y_size', 'min',
                'max', 'mean']  # , 'median', 'std_dev', 'percentile_1', 'percentile_99']

    def __init__(self, netcdf_path=None, netcdf_dataset=None,
                 max_bytes=500000000):
        '''
        DataStats Constructor
        Parameter:
            netcdf_path - string representing path to NetCDF file or URL for an OPeNDAP endpoint
            max_bytes - maximum number of bytes to pull into memory
        '''
        assert netcdf_dataset or netcdf_path, 'Either netcdf_dataset or netcdf_path must be defined'
        assert not (
                netcdf_dataset and netcdf_path), 'netcdf_dataset and netcdf_path cannot both be defined'

        netcdf_path = os.path.abspath(netcdf_path) if netcdf_path else None
        netcdf_dataset = netcdf_dataset or netCDF4.Dataset(netcdf_path, 'r')

        # Find variable with "grid_mapping" attribute - assumed to be 2D data
        # variable
        try:
            self.data_variable = [variable for variable in netcdf_dataset.variables.values(
            ) if hasattr(variable, 'grid_mapping')][0]
        except:
            raise Exception(
                'Unable to determine data variable (must have "grid_mapping" attribute')

        self._data_stats = {}
        self._data_stats['nc_path'] = netcdf_path or netcdf_dataset.filepath()
        self._data_stats['data_type'] = str(self.data_variable.dtype)
        self._data_stats['nodata_value'] = self.data_variable._FillValue

        shape = self.data_variable.shape
        # Array is ordered YX
        self._data_stats['x_size'] = shape[1]
        self._data_stats['y_size'] = shape[0]

        length_read = 0
        weighted_mean = 0.0

        for piece_array, _piece_offsets in array_pieces(
                self.data_variable, max_bytes=max_bytes):

            if isinstance(piece_array, np.ma.core.MaskedArray):
                piece_array = piece_array.data

            # Discard all no-data elements
            piece_array = np.array(
                piece_array[piece_array != self.data_variable._FillValue])

            piece_size = len(piece_array)

            if piece_size:
                try:
                    self._data_stats['min'] = min(
                        self._data_stats['min'], np.nanmin(piece_array))
                except:
                    self._data_stats['min'] = np.nanmin(piece_array)

                try:
                    self._data_stats['max'] = max(
                        self._data_stats['max'], np.nanmax(piece_array))
                except:
                    self._data_stats['max'] = np.nanmax(piece_array)

                weighted_mean += np.nanmean(piece_array) * piece_size
                length_read += piece_size
            # ==================================================================
            # else:
            #     print 'Empty array'
            # ==================================================================

        self._data_stats['mean'] = weighted_mean / length_read if length_read else None

        # ===================================================================
        # #TODO: Implement something clever for these
        # self._data_stats['median'] = np.NaN
        # self._data_stats['std_dev'] = np.NaN
        # self._data_stats['percentile_1'] = np.NaN
        # self._data_stats['percentile_99'] = np.NaN
        # ===================================================================

    # TODO: Do something nicer than this to get at the values, A property
    # might be good.
    def value(self, key):
        return self._data_stats[key]


def main():
    print(','.join(DataStats.key_list))
    for netcdf_path in sys.argv[1:]:
        datastats = DataStats(netcdf_path)
        print(','.join([str(datastats.value(key)) for key in DataStats.key_list]))


if __name__ == '__main__':
    main()
