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
'''
Created on 19 Oct. 2018

@author: Alex Ip - Geoscience Australia
'''
import sys
import argparse
import re
from datetime import datetime
from geophys_utils import CSWUtils
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

def date_string2datetime(date_string):
    '''
    Helper function to convert date string in one of several possible formats to a datetime object
    @param date_string: date string in one of several possible formats
    
    @return datetime object
    '''
    #?
    DATE_FORMAT_LIST = ['%Y%m{}', '%Y-%m-{}', '{}/%m/%y', '{}/%m/%Y']
    #If there is a date_string (start date or end date from argparse) then for each format type
    # in the DATE FORMAT LIST try to use the datetime.strptime method.
    #datetime.strptime() returns a datetime variable from the input parsed into the correct format.
    #  OR does it check if it is the correct format?
    datetime_result = None
    
    if date_string:
        for format_string in DATE_FORMAT_LIST:
            try:
                datetime_result = datetime.strptime(date_string, format_string)
                break
            except ValueError:
                pass
            
    #if successful return the input as a datetime class object or None.
    return datetime_result

def main():
    '''
    Main function to take command line parameters, perform CSW query and print required output
    '''
    def quote_delimitedtext(text, delimiter, quote_char='"'):
        '''
        Helper function to quote text containing delimiters or whitespace
        '''
        if delimiter in text or quote_char in text or re.search('\s', text):
            if delimiter == ',': # Use double quote to escape quote character for CSV
                return quote_char + text.replace(quote_char, quote_char + quote_char) + quote_char
            else: # Use backslash to escape quote character for tab or space delimited text
                return quote_char + text.replace(quote_char, '\\' + quote_char) + quote_char
        else:
            return text
            
        
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-k", "--keywords", help="comma-separated list of keywords for search", type=str)
    parser.add_argument("-t", "--titlewords", help="comma-separated list of titlewords for search", type=str)
    parser.add_argument("-a", "--anytext", help="comma-separated list of text snippets for search", type=str)
    parser.add_argument("-b", "--bounds", help='comma-separated <minx>,<miny>,<maxx>,<maxy> ordinates of bounding box for search. N.B: A leading "-" sign on this list should NOT be preceded by a space',
                        type=str)
    parser.add_argument("-c", "--crs", help='coordinate reference system for bounding box coordinates for search. Default determined by settings file.',
                        type=str)
    parser.add_argument("-s", "--start_date", help="start date for search", type=str)
    parser.add_argument("-e", "--end_date", help="end date for search", type=str)
    # Parameters to define output
    parser.add_argument("-p", "--protocols", help='comma-separated list of distribution protocols for output. Default determined by settings file, "*" = wildcard, "None" = no distributions.', type=str)
    parser.add_argument("-f", "--fields", help='comma-separated list of fields for output. Default determined by settings file, "*" = wildcard.', type=str)
    parser.add_argument("-d", "--delimiter", help='field delimiter for output. Defaults to "\t"', type=str)
    parser.add_argument("-u", "--urls", help="CSW URL(s) to query (comma separated list). Default determined by settings file.", type=str)
    parser.add_argument("-m", "--max_results", help="Maximum number of records to return", type=int)
    parser.add_argument('-r', '--header_row', action='store_const', const=True,
                        help='show header row. Default determined by settings file')
    parser.add_argument('-l', '--get_layers', action='store_const', const=True,
                        help='get WMS/WCS layer names. Default determined by settings file')
    parser.add_argument('--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("-y", "--types", help="comma-separated list of record types for search", type=str)
    
    args = parser.parse_args()


#CONVERTING INPUT TO CORRECT DATA TYPES AND FORMATS

    # Convert string to list of floats
    # if the user inputs an argument for bounds, then convert this string to a list of floats
    #  whereas each list object is split at the commas.
    # If the user does not call the bounds argument, then don't use it.
    if args.bounds:
        bounds = [float(ordinate) for ordinate in args.bounds.split(',')]
    else:
        bounds = None

    # Convert string to datetime
    start_date = date_string2datetime(args.start_date)

    # Convert string to datetime
    end_date = date_string2datetime(args.end_date)
    
    #creatse a CSW object and populates the parameters with argparse inputs

    url_list = ([url.strip() for url in args.urls.split(',')] if args.urls else None)

    cswu = CSWUtils(url_list,
                    debug=args.debug)

    # If there is a protocol list, then create a list of protocols that are split at the comma, use defaults if there isn't
    # Replace "None" with empty string
    protocol_list = (([protocol.strip().lower().replace('none', '') for protocol in args.protocols.split(',')] 
                     if args.protocols is not None else None) 
                     or cswu.settings['OUTPUT_DEFAULTS']['DEFAULT_PROTOCOLS'])
    
    # Allow wildcard - protocol_list=None means show all protocols
    if '*' in protocol_list:
        protocol_list = None

    # formatting the output for fields
    field_list = ([field.strip().lower() for field in args.fields.split(',')] if args.fields else None) or cswu.settings['OUTPUT_DEFAULTS']['DEFAULT_FIELDS']
    # Allow wildcard - field_list=None means show all fields
    if '*' in field_list:
        field_list = None
        
    # Set default delimiter to tab character
    delimiter = args.delimiter or cswu.settings['OUTPUT_DEFAULTS']['DEFAULT_DELIMITER']
    
    record_generator = cswu.query_csw(keyword_list=args.keywords,
                                 anytext_list=args.anytext,
                                 titleword_list=args.titlewords,
                                 bounding_box=bounds,
                                 start_datetime=start_date,
                                 stop_datetime=end_date,
                                 record_type_list=args.types,
                                 max_total_records=args.max_results,
                                 get_layers=args.get_layers
                                 )
    header_row_required = (cswu.settings['OUTPUT_DEFAULTS']['DEFAULT_SHOW_HEADER_ROW'] 
                           if args.header_row is None 
                           else args.header_row)

    # Print results
    header_printed = False
    distribution_count = 0
    for distribution in cswu.get_distributions(protocol_list, record_generator):
    #for distribution in cswu.get_netcdf_urls(record_generator):
        distribution_count += 1
        
        # Print header if required
        if header_row_required and not header_printed:
            print(delimiter.join([field
                                  for field in (field_list or sorted(distribution.keys()))
                                  ]
                                 )
                  )
            header_printed = True;
        
        # Decode and quote fields if required
        print(delimiter.join([quote_delimitedtext(str(distribution.get(field) or ''), delimiter)
                              for field in (field_list or sorted(distribution.keys()))
                              ]
                             )
              )

    logger.debug('{} results found.'.format(distribution_count))

if __name__ == '__main__':
    if not logger.handlers:
        # Set handler for root logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        #console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        
    main()