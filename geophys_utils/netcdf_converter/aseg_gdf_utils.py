'''
Functions to work with ASEG-GDF format string
Refer to https://www.aseg.org.au/sites/default/files/pdf/ASEG-GDF2-REV4.pdf for further information

Created on 19 Jun. 2018

@author: u76345
'''

import re
import numpy as np
from collections import OrderedDict
from math import ceil, log10
import logging
import warnings

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Logging level for this module

def dfrexp(f):
    '''
    Decimal version of frexp or np.frexp function to return mantissa & exponent
    @param f: Floating point scalar or array
    @return fman: Scalar or array decimal mantissa between 0.0 and 1.0 
    @return fexp: Scalar or array decimal exponent
    '''
    # Compute decimal exponent
    if type(f) == np.ndarray:
        fexp = np.zeros(shape=f.shape, dtype='int32')
        warnings.simplefilter("ignore") # Added to suppress warnings
        fexp[f != 0] = np.ceil(np.log10(np.abs(f[f != 0]))).astype('int32')
    else: # Scalar
        fexp = int(ceil(log10(abs(f)))) if f != 0 else 0
            
    # Compute decimal mantissa between 0.0 and 1.0
    fman = f/10.0**fexp
    
    logger.debug('fman: {}'.format(fman))
    logger.debug('fexp: {}'.format(fexp))
    
    return fman, fexp


# Approximate maximum number of significant decimal figures for each signed datatype
SIG_FIGS = OrderedDict([('uint8', 4), # 128
                        ('uint16', 10), # 32768
                        ('uint32', 19), # 2147483648 - should be 9, but made 10 because int64 is unsupported
                        ('uint64', 30), # 9223372036854775808 - Not supported in netCDF3 or netCDF4-Classic
                        ('int8', 2), # 128
                        ('int16', 4), # 32768
                        ('int32', 10), # 2147483648 - should be 9, but made 10 because int64 is unsupported
                        ('int64', 19), # 9223372036854775808 - Not supported in netCDF3 or netCDF4-Classic
                        # https://en.wikipedia.org/wiki/Floating-point_arithmetic#IEEE_754:_floating_point_in_modern_computers
                        ('float32', 7), # 7.2
                        ('float64', 35) # 15.9 - should be 16, but made 35 to support unrealistic precision specifications
                        ]
                       )

DTYPE_REDUCTION_LISTS = [['int64', 'int32', 'int16', 'int8'], # Integer dtypes
                         ['uint64', 'uint32', 'uint16', 'uint8'], # Unsigned integer dtypes
                         ['float64', 'float32'] #, 'int16', 'int8'] # Floating point dtypes - do NOT try casting to integer types
                         ]
    
ASEG_DTYPE_CODE_MAPPING = {'uint8': 'I',
                           'uint16': 'I',
                           'uint32': 'I',
                           'uint64': 'I',
                           'int8': 'I',
                           'int16': 'I',
                           'int32': 'I',
                           'int64': 'I',
                           'float32': 'E', # real in exponent form
                           'float64': 'D', # double precision real in exponent form
                           'str': 'A'
                           }

def decode_aseg_gdf_format(aseg_gdf_format):
    '''
    Function to decode ASEG-GDF format string
    @param aseg_gdf_format: ASEG-GDF format string

    @return columns: Number of columns (i.e. 1 for 1D data, or read from format string for 2D data)
    @return aseg_dtype_code: ASEG-GDF data type character, e.g. "F" or "I"
    @return width_specifier:  Width of field in number of characters read from format string
    @return decimal_places: Number of fractional digits read from format string    
    '''
    if not aseg_gdf_format:
        raise BaseException('No ASEG-GDF format string to decode')  

    match = re.match('(\d+)*(\w)(\d+)\.*(\d+)*', aseg_gdf_format)
    
    if not match:
        raise BaseException('Invalid ASEG-GDF format string {}'.format(aseg_gdf_format))  
      
    columns = int(match.group(1)) if match.group(1) is not None else 1
    aseg_dtype_code = match.group(2).upper()
    width_specifier = int(match.group(3))
    decimal_places = int(match.group(4)) if match.group(4) is not None else 0
    
    logger.debug('aseg_gdf_format: {}, columns: {}, aseg_dtype_code: {}, width_specifier: {}, decimal_places: {}'.format(aseg_gdf_format, 
                                                                                                                      columns, 
                                                                                                                      aseg_dtype_code, 
                                                                                                                      width_specifier, 
                                                                                                                      decimal_places
                                                                                                                      )
                 ) 
    return columns, aseg_dtype_code, width_specifier, decimal_places  

def aseg_gdf_format2dtype(aseg_gdf_format):
    '''
    Function to return Python data type string and other precision information from ASEG-GDF format string
    @param aseg_gdf_format: ASEG-GDF format string

    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or read from format string for 2D data)
    @return width_specifier:  Width of field in number of characters read from format string
    @return decimal_places: Number of fractional digits read from format string    
    '''
    columns, aseg_dtype_code, width_specifier, decimal_places = decode_aseg_gdf_format(aseg_gdf_format)
    dtype = None # Initially unknown datatype
    
    # Determine type and size for required significant figures
    # Integer type - N.B: Only signed types available
    if aseg_dtype_code == 'I':
        assert not decimal_places, 'Integer format cannot be defined with fractional digits'
        for test_dtype, sig_figs in SIG_FIGS.items():
            if test_dtype.startswith('int') and sig_figs >= width_specifier:
                dtype = test_dtype
                break
        assert dtype, 'Invalid width_specifier of {}'.format(width_specifier)     
    
    # Floating point type - use approximate sig. figs. to determine length
    #TODO: Remove 'A' after string field handling has been properly implemented
    elif aseg_dtype_code in ['D', 'E', 'F']: # Floating point
        for test_dtype, sig_figs in SIG_FIGS.items():
            if test_dtype.startswith('float') and sig_figs >= width_specifier-2: # Allow for sign and decimal place
                dtype = test_dtype
                break
        assert dtype, 'Invalid floating point format of {}.{}'.format(width_specifier, decimal_places)                                    
    
    elif aseg_dtype_code == 'A':
        assert not decimal_places, 'String format cannot be defined with fractional digits'
        dtype = '<U{}'.format(width_specifier) # Numpy fixed-length string type
        
    else:
        raise BaseException('Unhandled ASEG-GDF dtype code {}'.format(aseg_dtype_code))
    
    logger.debug('aseg_dtype_code: {}, columns: {}, width_specifier: {}, decimal_places: {}'.format(dtype, 
                                                                                                 columns, 
                                                                                                 width_specifier, 
                                                                                                 decimal_places
                                                                                                 )
                 ) 
    return dtype, columns, width_specifier, decimal_places


def variable2aseg_gdf_format(array_variable, decimal_places=None):
    '''
    Function to return ASEG-GDF format string and other info from data array or netCDF array variable
    @param array_variable: data array or netCDF array variable
    @param decimal_places: Number of decimal places to respect, or None for value derived from datatype and values
    
    @return aseg_gdf_format: ASEG-GDF format string
    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or second dimension size for 2D data)
    @return width_specifier: Width of field in number of characters
    @return decimal_places: Number of fractional digits (derived from datatype sig. figs - width_specifier)
    @param python_format: Python Formatter string for fixed-width output
    '''
    if len(array_variable.shape) <= 1: # 1D variable or scalar
        columns = 1
    elif len(array_variable.shape) == 2: # 2D variable
        columns = array_variable.shape[1]
    else:
        raise BaseException('Unable to handle arrays with dimensionality > 2')
        
    data_array = array_variable[:]
    
    # Try to determine the dtype string from the variable and data_array type
    if not len(array_variable.shape): # Scalar
        dtype = type(data_array).__name__
        if dtype == 'str':
            width_specifier = len(data_array) + 1
            decimal_places = 0 
        elif dtype == 'ndarray': # Single-element array
            dtype = str(array_variable.dtype)
            data = np.asscalar(data_array)

            sig_figs = SIG_FIGS[dtype] + 1 # Look up approximate significant figures and add 1
            sign_width = 1 if data < 0 else 0
            integer_digits = ceil(log10(np.abs(data) + 1.0))
    else: # Array
        dtype = str(array_variable.dtype)
        if dtype in ['str', "<class 'str'>"]: # String array or string array variable
            dtype = 'str'
            width_specifier = max([len(string.strip()) for string in data_array]) + 1
            decimal_places = 0
            
        else:  # Numeric datatype array
            # Include fill value if required
            if type(data_array) == np.ma.core.MaskedArray:
                logger.debug('Array is masked. Including fill value.')
                data_array = data_array.data
                
            sig_figs = SIG_FIGS[dtype] + 1 # Look up approximate significant figures and add 1
            sign_width = 1 if np.nanmin(data_array) < 0 else 0
            integer_digits = ceil(log10(np.nanmax(np.abs(data_array)) + 1.0))
    
    aseg_dtype_code = ASEG_DTYPE_CODE_MAPPING.get(dtype)
    assert aseg_dtype_code, 'Unhandled dtype {}'.format(dtype)
    
    if aseg_dtype_code == 'I': # Integer
        decimal_places = 0
        width_specifier = integer_digits + sign_width + 1
        aseg_gdf_format = 'I{}'.format(width_specifier)
        python_format = '{' + ':>{:d}.{:d}f'.format(width_specifier, decimal_places) + '}'

    elif aseg_dtype_code in ['F', 'D', 'E']: # Floating point
        # If array_variable is a netCDF variable with a "format" attribute, use stored format string to determine decimal_places
        if decimal_places is not None:
            decimal_places = min(decimal_places, abs(sig_figs-integer_digits))
            logger.debug('decimal_places set to {} from decimal_places {}'.format(decimal_places, decimal_places))
        elif hasattr(array_variable, 'aseg_gdf_format'): 
            _columns, _aseg_dtype_code, _integer_digits, decimal_places = decode_aseg_gdf_format(array_variable.aseg_gdf_format)
            decimal_places = min(decimal_places, abs(sig_figs-integer_digits))
            logger.debug('decimal_places set to {} from variable attribute aseg_gdf_format {}'.format(decimal_places, array_variable.aseg_gdf_format))
        else: # No aseg_gdf_format variable attribute
            decimal_places = abs(sig_figs-integer_digits) # Allow for full precision of datatype
            logger.debug('decimal_places set to {} from sig_figs {} and integer_digits {}'.format(decimal_places, sig_figs, integer_digits))
        
        width_specifier = min(sign_width + integer_digits + decimal_places + 2,
                              sign_width + sig_figs + 2
                              )
                                  
        aseg_gdf_format = '{}{}.{}'.format(aseg_dtype_code, width_specifier, decimal_places)
        if aseg_dtype_code == 'F': # Floating point notation
            python_format = '{' + ':>{:d}.{:d}f'.format(width_specifier, decimal_places) + '}' # Add 1 to width for decimal point
        else: # Exponential notation for 'D' or 'E'
            python_format = '{' + ':>{:d}.{:d}E'.format(width_specifier, decimal_places) + '}' # Add 1 to width for decimal point

    elif aseg_dtype_code == 'A': # String
        if hasattr(array_variable, 'aseg_gdf_format'):
            _columns, _aseg_dtype_code, width_specifier, decimal_places = decode_aseg_gdf_format(array_variable.aseg_gdf_format)
            aseg_gdf_format = array_variable.aseg_gdf_format
        else:
            aseg_gdf_format = 'A{}'.format(width_specifier)
            
        python_format = '{' + ':>{:d}s'.format(width_specifier) + '}'
    else:
        raise BaseException('Unhandled ASEG-GDF dtype code {}'.format(aseg_dtype_code))
    
    # Pre-pend column count to start of aseg_gdf_format
    if columns > 1:
        aseg_gdf_format = '{}{}'.format(columns, aseg_gdf_format)
        
    return aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format


def fix_field_precision(data_array, current_dtype, decimal_places, no_data_mask=[], fill_value=None):
    '''
    Function to return revised ASEG-GDF format string and other info from data array or netCDF array variable
    after correcting datatype for excessive precision specification, or None if there is no precision change.
    Arrays are copied to smaller representations and then the difference with the original is checked to
    ensure that any difference is less than precision of the specified number of fractional digits.
    Note that fill_value is also considered but potentially modified only if data precision is changed
    @param data_array: data array - assumed to be of dtype float64 for raw data
    @param current_dtype: Current data type string, e.g. int8 or float32
    @param decimal_places: Number of fractional digits for precision checking
    @param fill_value: fill value or None
    
    Returns None if no precision change required.
    @return aseg_gdf_format: ASEG-GDF format string
    @return dtype: Data type string, e.g. int8 or float32
    @return columns: Number of columns (i.e. 1 for 1D data, or second dimension size for 2D data)
    @return width_specifier:  Width of field in number of characters
    @return decimal_places: Number of fractional digits (derived from datatype sig. figs - width_specifier)
    @return python_format: Python Formatter string for fixed-width output
    @return fill_value: Potentially modified fill value
    '''
    logger.debug('data_array: {}, current_dtype: {}, decimal_places: {}'.format(data_array, current_dtype, decimal_places))
    
    try:
        data_mantissa, data_exponent = dfrexp(data_array)
    except:
        logger.debug('Unable to compute data_mantissa & data_exponent')
        return
    
    for dtype_reduction_list in DTYPE_REDUCTION_LISTS:
        try:
            current_dtype_index = dtype_reduction_list.index(current_dtype)

            # Try types from smallest to largest
            for smaller_dtype in dtype_reduction_list[:current_dtype_index:-1]:                
                smaller_array = data_array.astype(smaller_dtype)
                difference_array = data_array - smaller_array
                logger.debug('current_dtype: {}\nsmaller_dtype: {}\narray_variable\n{}\nsmaller_array\n{}\n\
difference_array\n{}\ndecimal_places: {}\ndifference count: {}\ndifference values: '.format(current_dtype, 
                                                                                               smaller_dtype, 
                                                                                               data_array, 
                                                                                               smaller_array, 
                                                                                               difference_array, 
                                                                                               decimal_places, 
                                                                                               np.count_nonzero(difference_array >= pow(10, -decimal_places)), 
                                                                                               difference_array[difference_array != 0]
                                                                                               )
                      )
                logger.debug('Maximum error converting from {} to {}: {}'.format(current_dtype,
                                                                                 smaller_dtype,
                                                                                 np.nanmax(np.abs(difference_array))))
                smaller_mantissa, smaller_exponent = dfrexp(smaller_array)
                if np.any(np.logical_or((smaller_exponent != data_exponent), 
                                        (np.abs(data_mantissa - smaller_mantissa) >= pow(10, -decimal_places)))):
                    # Differences found - try larger datatype
                    continue
                else:
                    logger.debug('Maximum mantissa difference: {}'.format(np.nanmax(np.abs(data_mantissa - smaller_mantissa))))
                    logger.debug('Maximum exponent difference: {}'.format(np.nanmax(np.abs(smaller_exponent - data_exponent))))
                    aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format = variable2aseg_gdf_format(smaller_array, decimal_places)

                    if fill_value is not None:
                        # Use reduced precision fill_value if available and unambiguous
                        if np.any(no_data_mask):
                            reduced_precision_fill_value = data_array[no_data_mask][0] 
                            
                            # Check for ambiguity introduced by reduced precision
                            if np.any(data_array[~no_data_mask] == reduced_precision_fill_value):
                                logger.debug('Reduced precision fill value of {} is ambiguous'.format(reduced_precision_fill_value))
                                continue # Can't use this datatype - try the next bigger one
                            elif fill_value != reduced_precision_fill_value:
                                logger.debug('fill_value precision reduced from {} to {}'.format(fill_value, 
                                                                                                reduced_precision_fill_value)
                                                                                                )
                                fill_value = reduced_precision_fill_value
                        
                        fill_value = truncate(fill_value, data_array, no_data_mask, width_specifier, decimal_places)
                        
                    return aseg_gdf_format, dtype, columns, width_specifier, decimal_places, python_format, fill_value

            
        except ValueError: # current_dtype not in dtype_reduction_list
            continue


def truncate(fill_value, data_array, no_data_mask, width_specifier, decimal_places):
    '''
    Function to truncate fill_value to <width_specifier>.<decimal_places> rather than rounding for neater output later on

    @param fill_value: Original fill value
    @param data_array: Array containing data
    @param no_data_mask: Boolean mask array (true for valid data)
    @param width_specifier:  Width of field in number of characters
    @param decimal_places: Number of fractional digits for precision checking
    
    @return fill_value: Potentially modified fill value
    '''
    try:
        truncated_fill_value = None
        
        integer_digits = width_specifier - decimal_places
        if fill_value < 0:
            integer_digits -= 1 # Allow for sign
        if decimal_places > 0:
            integer_digits -= 1 # Allow for decimal point
            
        fill_value_str = (fill_value)
        assert 'e' not in fill_value_str.lower(), 'Unable to truncate value in exponential notation'
        pattern = re.compile('(-?)\d*?(\d{0,' + '{}'.format(integer_digits) + '}\.\d{0,' + '{}'.format(decimal_places) + '})')
        search = re.search(pattern, fill_value_str)
        truncated_fill_value = float(search.group(1)+search.group(2))
        # Check for any ambiguity introduced by truncation
        assert not np.any(data_array[~no_data_mask] == truncated_fill_value), 'Truncated fill value of {} is ambiguous'.format(truncated_fill_value)
        if fill_value != truncated_fill_value:
            logger.debug('fill_value truncated from {} to {}'.format(fill_value, truncated_fill_value))
        return truncated_fill_value
    except Exception as e:
        logger.debug('Unable to truncate fill value from {} to {} ({}). Keeping original value.'.format(fill_value, truncated_fill_value, e))
        return fill_value
