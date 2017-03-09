"""
Unit tests for geophys_utils._array_pieces module

Created on 15/11/2016

@author: Alex Ip
"""
import unittest
import numpy as np
from geophys_utils._array_pieces import array_pieces

class TestArrayPieces(unittest.TestCase):
    """Unit tests for geophys_utils._array_pieces module."""
    
    def test_array_pieces(self):
        print 'Testing array_pieces function'
        # create 100 x 100 array with unique elements
        test_array = np.reshape(np.arange(0, 10000, dtype='int16'), (100,100))
        sixteenth_bytes = test_array.dtype.itemsize * reduce(lambda x, y: x * y / 16, test_array.shape)
        overlap = 10
        
        print'\tTesting single piece'
        array_pieces_results = {array_offset: piece_array for piece_array, array_offset in array_pieces(test_array)}

        assert len(array_pieces_results) == 1, 'Whole array not returned for large max_bytes'
        assert not np.any(test_array - array_pieces_results.values()[0]), 'Array contents changed'
        
        print'\tTesting sixteenth arrays'
        array_pieces_results = {array_offset: piece_array for piece_array, array_offset in array_pieces(test_array,
                                                                                                        max_bytes=sixteenth_bytes)}
        assert len(array_pieces_results) == 16, '16 arrays not returned for max_bytes=%d' % sixteenth_bytes
        for array_offset in sorted(array_pieces_results.keys()):
            piece_array = array_pieces_results[array_offset]
            
            assert piece_array.dtype.itemsize * reduce(lambda x, y: x * y / 16, piece_array.shape) <= sixteenth_bytes, 'piece_array too large'
            
            expected_shape = tuple([test_array.shape[dim_index] / 4 
                                    for dim_index in range(2)])
            assert piece_array.shape == expected_shape, 'piece_array is wrong shape'

            slices = [slice(array_offset[dim_index], 
                            array_offset[dim_index]+piece_array.shape[dim_index]
                            ) for dim_index in range(2)]
            assert not np.any(piece_array - test_array[slices]), 'Array contents changed for array piece at %s' % array_offset
        
        print'\tTesting sixteenth arrays with overlap'
        array_pieces_results = {array_offset: piece_array for piece_array, array_offset in array_pieces(test_array,
                                                                                                        max_bytes=sixteenth_bytes,
                                                                                                        overlap=overlap)}

        assert len(array_pieces_results) == 16, '16 arrays not returned for overlap=%d, max_bytes=%d' % (overlap,
                                                                                                                  sixteenth_bytes)
        for array_offset in sorted(array_pieces_results.keys()):
            piece_array = array_pieces_results[array_offset] 

            expected_shape = tuple([test_array.shape[dim_index] / 4 + overlap 
                                    if array_offset[dim_index] == 0 
                                    or test_array.shape[dim_index] - array_offset[dim_index] - overlap <= test_array.shape[dim_index] / 4
                                    else test_array.shape[dim_index] / 4 + 2 * overlap 
                                    for dim_index in range(2)]
                                    )
            assert piece_array.shape == expected_shape, 'piece_array is wrong shape'

            slices = [slice(array_offset[dim_index], 
                            array_offset[dim_index]+piece_array.shape[dim_index]
                            ) for dim_index in range(2)]
            assert not np.any(piece_array - test_array[slices]), 'Array contents changed for array piece at %s' % array_offset
        


# Define test suites
def test_suite():
    """Returns a test suite of all the tests in this module."""

    test_classes = [TestArrayPieces]

    suite_list = map(unittest.defaultTestLoader.loadTestsFromTestCase,
                     test_classes)

    suite = unittest.TestSuite(suite_list)

    return suite


# Define main function
def main():
    unittest.TextTestRunner(verbosity=2).run(test_suite())

if __name__ == '__main__':
    main()
