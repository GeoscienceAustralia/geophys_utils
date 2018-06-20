'''
Functions to work with ASEG-GDF format string

Created on 19 Jun. 2018

@author: u76345
'''

import re
import numpy as np
from collections import OrderedDict
from math import log10, ceil

# Approximate maximum number of significant decimal figures for each signed datatype
SIG_FIGS = OrderedDict([('int8', 2), # 128
                        ('int16', 4), # 32768
                        ('int32', 10), # 2147483648 - should be 9, but made 10 because int64 is unsupported
                        ('int64', 19), # 9223372036854775808
                        # https://en.wikipedia.org/wiki/Floating-point_arithmetic#IEEE_754:_floating_point_in_modern_computers
                        ('float32', 7), # 7.2
                        ('float64', 20) # 15.9 - should be 16, but made 20 to support unrealistic precision specifications
                        ]
                       )

def aseg_gdf_format2dtype(aseg_gdf_format):
    '''
    Function to return data type string from ASEG-GDF format string
    Refer to https://www.aseg.org.au/sites/default/files/pdf/ASEG-GDF2-REV4.pdf for further information
    @param aseg_gdf_format: ASEG-GDF format string

    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or read from format string for 2D data)
    @return integer_digits: Number of integer digits read from format string
    @return fractional_digits: Number of fractional digits read from format string
    '''
    if not aseg_gdf_format:
        raise BaseException('No ASEG-GDF format string to decode')  

    match = re.match('(\d+)*(\w)(\d+)\.*(\d+)*', aseg_gdf_format)
    
    if not match:
        raise BaseException('Invalid ASEG-GDF format string {}'.format(aseg_gdf_format))  
      
    columns = match.group(1) or 1
    dtype_code = match.group(2).upper()
    integer_digits = int(match.group(3))
    fractional_digits = int(match.group(4)) if match.group(4) is not None else 0
    dtype = None # Initially unknown datatype
    
    # Determine type and size for required significant figures
    # Integer type - N.B: Only signed types available
    if dtype_code == 'I':
        assert not fractional_digits, 'Integer format cannot be defined with fractional digits'
        for test_dtype, sig_figs in SIG_FIGS.items():
            if test_dtype.startswith('int') and sig_figs >= integer_digits:
                dtype = test_dtype
                break
        assert dtype, 'Invalid integer length of {}'.format(integer_digits)     
    
    # Floating point type - use approximate sig. figs. to determine length
    elif dtype_code in ['D', 'E', 'F']: 
        for test_dtype, sig_figs in SIG_FIGS.items():
            if test_dtype.startswith('float') and sig_figs >= integer_digits + fractional_digits:
                dtype = test_dtype
                break
        assert dtype, 'Invalid floating point format of {}.{}'.format(integer_digits, fractional_digits)                                    
    
    return dtype, columns, integer_digits, fractional_digits


def dtype2aseg_gdf_format(array_variable):
    '''
    Function to return ASEG-GDF format string and other info from data array or netCDF array variable
    Refer to https://www.aseg.org.au/sites/default/files/pdf/ASEG-GDF2-REV4.pdf for further information
    @param array_variable: data array or netCDF array variable
    
    @return aseg_gdf_format: ASEG-GDF format string
    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or second dimension size for 2D data)
    @return integer_digits: Number of integer digits (derived from maximum value)
    @return fractional_digits: Number of fractional digits (derived from datatype sig. figs - integer_digits)
    @param python_format: Python Formatter string for fixed-width output
    '''
    if len(array_variable.shape) == 1: # 1D variable
        columns = 1
    elif len(array_variable.shape) == 2: # 2D variable
        columns = array_variable.shape[1]
    else:
        raise BaseException('Unable to handle arrays with dimensionality > 3')
        
    dtype = str(array_variable.dtype)
            
    sig_figs = SIG_FIGS[dtype] + 1 # Look up approximate significant figures and add 1
    sign_width = 1 if np.nanmin(array_variable[:]) < 0 else 0
    integer_digits = ceil(log10(np.nanmax(np.abs(array_variable[:])) + 1.0))
    
    if dtype.startswith('int'):
        fractional_digits = 0
        aseg_gdf_format = 'I{}'.format(integer_digits)
        python_format = '{' + ':{:d}.{:d}f'.format(sign_width+integer_digits, fractional_digits) + '}'
    elif dtype.startswith('float'):
        # If array_variable is a netCDF variable with a "format" attribute, use stored format string
        if hasattr(array_variable, 'aseg_gdf_format'): 
            _dtype, _columns, _integer_digits, fractional_digits = aseg_gdf_format2dtype(array_variable.aseg_gdf_format)
            fractional_digits = min(fractional_digits, sig_figs-integer_digits)
        else: # No fra
            fractional_digits = sig_figs - integer_digits
            
        aseg_gdf_format = 'F{}.{}'.format(integer_digits, fractional_digits)
        python_format = '{' + ':{:d}.{:d}f'.format(sign_width+sig_figs+1, fractional_digits) + '}' # Add 1 to width for decimal point
    
    # Pre-pend column count to start of aseg_gdf_format
    if columns > 1:
        aseg_gdf_format = '{}{}'.format(columns, aseg_gdf_format)
        
    return aseg_gdf_format, dtype, columns, integer_digits, fractional_digits, python_format


def fix_field_precision(array_variable, current_dtype, fractional_digits):
    '''
    Function to return revised ASEG-GDF format string and other info from data array or netCDF array variable
    after correcting excessive precision specification, or None if there is no precision change.
    @param array_variable: data array or netCDF array variable - assumed to be of dtype float64 for raw data
    @param current_dtype: Current data type string, e.g. int8 or float32
    @param fractional_digits: Number of fractional digits for precision checking
    
    Returns None if no precision change required.
    @return aseg_gdf_format: ASEG-GDF format string
    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or second dimension size for 2D data)
    @return integer_digits: Number of integer digits (derived from maximum value)
    @return fractional_digits: Number of fractional digits (derived from datatype sig. figs - integer_digits)
    @param python_format: Python Formatter string for fixed-width output
    '''
    dtype_reduction_lists = [['int64', 'int32', 'int16', 'int8'], # Integer dtypes
                             ['float64', 'float32'] # Floating point dtypes
                             ]
    
    for dtype_reduction_list in dtype_reduction_lists:
        try:
            current_dtype_index = dtype_reduction_list.index(current_dtype)
            # Try types from smallest to largest
            for smaller_dtype in dtype_reduction_list[:current_dtype_index:-1]:                
                smaller_array = array_variable[:].astype(smaller_dtype)
                difference_array = array_variable[:] - smaller_array
                #===============================================================
                # print('current_dtype:', current_dtype,
                #       '\nsmaller_dtype:', smaller_dtype,
                #       '\narray_variable\n', 
                #       array_variable[:], 
                #       '\nsmaller_array\n', 
                #       smaller_array,
                #       '\ndifference_array\n', 
                #       difference_array,
                #       '\nfractional_digits:', fractional_digits,
                #       '\ndifference count:', np.count_nonzero(difference_array >= pow(10, -fractional_digits)),
                #       '\ndifference values: ', difference_array[difference_array != 0]
                #       )
                #===============================================================
                if np.count_nonzero(np.abs(difference_array) >= pow(10, -fractional_digits)):
                    # Differences found - try larger datatype
                    continue
                else:
                    reduced_precision_result = dtype2aseg_gdf_format(smaller_array)
                    # aseg_gdf_format, dtype, columns, integer_digits, fractional_digits, python_format
                    # Use the minimum number of decimal places: either specified or derived
                    reduced_precision_result[4] = min(reduced_precision_result[4], fractional_digits)
                    return reduced_precision_result

            
        except ValueError: # current_dtype_index not in dtype_reduction_list
            continue
   
