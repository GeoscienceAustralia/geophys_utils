'''
process_array.py
Created on 1 Sep,2016

@author: Alex Ip
'''
import sys
import netCDF4
import math
import itertools
from functools import reduce

def array_pieces(ndarray, overlap=0, max_bytes=None):
    '''
    Generator to return a series of numpy arrays less than max_bytes in size and the offset within the complete data from a NetCDF variable
    Parameters:
        ndarray: Numpy array or NetCDF array variable
        overlap: number of pixels to add to each edge
        max_bytes: Maximum number of bytes to retrieve. Defaults to 500,000,000 for NCI's OPeNDAP

    Yields:
        piece_array: array subset less than max_bytes in size
        array_offset: start indices of subset in whole array
    '''
    max_bytes = max_bytes or 500000000  # Defaults to 500,000,000 for NCI's OPeNDAP

    array_shape = ndarray.shape
    array_dimensions = len(array_shape)

    # Determine overall array size in bytes
    array_bytes = ndarray.dtype.itemsize * \
        reduce(lambda x, y: x * y, array_shape)

    if array_bytes > max_bytes:  # Multiple pieces required
        # Determine number of divisions in each axis required to keep pieces
        # under max_bytes in size
        axis_divisions = int(math.ceil(
            math.pow(math.ceil(array_bytes / float(max_bytes)), 1.0 / array_dimensions)))

        # Determine chunk size for pieces or default to natural divisions if no
        # chunking set
        try:
            chunking = ndarray.chunking() or (1, 1)
        except:  # Numpy arrays don't have chunking
            chunking = (1, 1)

        # Determine piece shape rounded down to chunk sizes
        piece_shape = [array_shape[index] / axis_divisions / chunking[index]
                       * chunking[index] for index in range(array_dimensions)]

        # Determine number of pieces in each axis - all elements should all be
        # the same, so this is probably redundant
        axis_pieces = [array_shape[index] / piece_shape[index]
                       for index in range(array_dimensions)]

        # Iterate over every piece of array
        for piece_indices in itertools.product(*[range(axis_pieces[dimension_index]) 
                                                 for dimension_index in range(array_dimensions)]):
            
            # Compute base start indices with no overlap
            start_indices = [piece_indices[dimension_index] * piece_shape[dimension_index]
                             for dimension_index in range(array_dimensions)]
            
            # Compute end indices plus overlap from start indices
            end_indices = [min(start_indices[dimension_index] + piece_shape[dimension_index] + overlap,
                               array_shape[dimension_index]) 
                           for dimension_index in range(array_dimensions)]
            
            # Subtract overlap from base start indices 
            start_indices = [max(0, start_indices[dimension_index] - overlap)
                             for dimension_index in range(array_dimensions)]
            
            array_slices = [slice(start_indices[dimension_index],
                                  end_indices[dimension_index])
                            for dimension_index in range(array_dimensions)]

            piece_array = ndarray[array_slices]
            yield piece_array, tuple(start_indices)

    else:  # Only one piece required
        yield ndarray[...], (0, 0)


def main():
    '''
    Main function for testing
    '''
    netcdf_path = sys.argv[1]
    netcdf_dataset = netCDF4.Dataset(netcdf_path)

    # Find variable with "grid_mapping" attribute - assumed to be 2D data
    # variable
    try:
        data_variable = [variable for variable in netcdf_dataset.variables.values(
        ) if hasattr(variable, 'grid_mapping')][0]
    except:
        raise Exception(
            'Unable to determine data variable (must have "grid_mapping" attribute')

    piece_count = 0
    for piece_array, array_offset in array_pieces(data_variable, overlap=0):
        piece_count += 1
        piece_bytes = data_variable.dtype.itemsize * \
            reduce(lambda x, y: x * y, piece_array.shape)
        print 'piece_array.shape = %s, array_offset = %s, piece_bytes = %d' % (piece_array.shape, array_offset, piece_bytes)

    print 'piece_count = %s' % piece_count

if __name__ == '__main__':
    main()
