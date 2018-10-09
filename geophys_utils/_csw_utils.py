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
Created on 23Feb.,2017

@author: u76345
'''
import re
import os
import copy
from datetime import datetime
from owslib import fes
import argparse
from owslib.csw import CatalogueServiceWeb
from owslib.wms import WebMapService
from owslib.wcs import WebCoverageService
import netCDF4
import yaml
from pprint import pformat
import logging
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

class CSWUtils(object):
    '''
    CSW query utilities
    '''

    def __init__(self, 
                 csw_url_list=None, 
                 timeout=None,
                 debug=False,
                 settings_path=None
                 ):
        '''
        Constructor for CSWUtils class
        @param csw_url_list: List of URLs for CSW services. Defaults to value of self.settings['DEFAULT_CSW_URLS']
        @param timeout: Timeout in seconds. Defaults to value of self.settings['DEFAULT_TIMEOUT']
        @param settings_path: Path to settings file defaults to csw_utils_settings.yml in package directory
        '''
        self._debug = None # Initialise private variable
        self.debug = debug # Set debug property
        
        package_dir = os.path.dirname(os.path.abspath(__file__))
        settings_path = settings_path or os.path.join(package_dir, 'csw_utils_settings.yml')
        self.settings = yaml.safe_load(open(settings_path))
        
        logger.debug('Settings: {}'.format(pformat(self.settings)))
        
        self.csw_url_list = csw_url_list or self.settings['DEFAULT_CSW_URLS']
        if type(self.csw_url_list) == str:
            self.csw_url_list = [self.csw_url_list]
            
        timeout = timeout or self.settings['DEFAULT_TIMEOUT']
        
        for key, value in iter(self.settings['ENVIRONMENT_VARIABLES'].items()):
            if value:
                os.environ[key] = value

        self.csw_list = []
        for csw_url in self.csw_url_list:
            try:
                self.csw_list.append(CatalogueServiceWeb(csw_url, timeout=timeout))
            except Exception as e:
                logger.warning('Unable to open CSW URL {}: {}'.format(csw_url, e))
                
        assert self.csw_list, 'No valid CSW URLs provided'

    def list_from_comma_separated_string(self, comma_separated_string):
        '''
        Helper function to return list of strings from a comma-separated string
        @param comma_separated_string: comma-separated string
        @return: list of strings
        '''
        return [keyword.strip() for keyword in comma_separated_string.split(',')]

    # A helper function for date range filtering
    def get_date_filter(self, start_datetime=None, stop_datetime=None, constraint='overlaps'):
        '''
        Helper function to return a list containing a pair of FES filters for a date range
        @param  start_datetime: datetime object for start of time period to search
        @param stop_datetime: datetime object for end of time period to search
        @param constraint: string value of either 'overlaps' or 'within' to indicate type of temporal search

        @return: list containing a pair of FES filters for a date range
        '''
        if start_datetime:
            start_date_string = start_datetime.date().isoformat()
        else:
            start_date_string = '1900-01-01' # Distant past

        if stop_datetime:
            stop_date_string = stop_datetime.date().isoformat()
        else:
            stop_date_string = '2100-01-01' # Distant future
            
        if constraint == 'overlaps':
            start_filter = fes.PropertyIsLessThanOrEqualTo(propertyname='ows100:TempExtent_begin', literal=stop_date_string)
            stop_filter = fes.PropertyIsGreaterThanOrEqualTo(propertyname='ows100:TempExtent_end', literal=start_date_string)
        elif constraint == 'within':
            start_filter = fes.PropertyIsGreaterThanOrEqualTo(propertyname='ows100:TempExtent_begin', literal=start_date_string)
            stop_filter = fes.PropertyIsLessThanOrEqualTo(propertyname='ows100:TempExtent_end', literal=stop_date_string)

        return [start_filter, stop_filter]
    
    def any_case_match_filters(self, propertyname, word_list):
        '''
        Helper function to return and/or filter from a word list for partial case insensitivity
        Work-around resolves issues where "matchCase=False" is ignored and keywords with whitespace characters are not found
        @param propertyname: String denoting CSW property to search
        @param word_list: List of strings for (partially) case-insensitive "and" search
        
        @return: List of FES filters expressing "and" query
        '''
        filter_list = []
        for word in word_list:
            word_variant_set = set([word, word.upper(), word.lower(), word.title()])
            
            if len(word_variant_set) == 1: # Use single filter if no "or" required
                filter_list.append(fes.PropertyIsLike(propertyname=propertyname, literal=word_variant_set.pop(), matchCase=False))
            else:
                filter_list.append(fes.Or([fes.PropertyIsLike(propertyname=propertyname, 
                                                              literal=re.sub('(\s)', '_', word_variant), 
                                                              matchCase=False
                                                              )
                                           for word_variant in set([word, word.upper(), word.lower(), word.title()])
                                           ]
                                          )
                                   )
                
        return filter_list

    def get_csw_records(self, 
                        fes_filters, 
                        max_query_records=None, 
                        max_total_records=None,
                        get_layers=False
                        ):
        '''
        Generator yeilding nested dicts containing information about each CSW query result record including distributions
        @param fes_filters: List of fes filters to apply to CSW query
        @param max_query_records: Maximum number of records to return per CSW query. Defaults to value of self.settings['DEFAULT_MAXRECORDS']
        @param max_total_records: Maximum total number of records to return. Defaults to value of self.settings['DEFAULT_MAXTOTALRECORDS']
        @param get_layers: Boolean flag indicating whether to get WMS/WCS layer names. Defaults to False
        '''
        max_query_records = max_query_records or self.settings['DEFAULT_MAXQUERYRECORDS']
        max_total_records = max_total_records or self.settings['DEFAULT_MAXTOTALRECORDS']
        
        uuid_list = []
        for csw in self.csw_list:
            logger.debug('Querying {}'.format(csw.url))

            startposition = 1 # N.B: This is 1-based, not 0-based
            record_count = 0
    
            # Keep querying until all results have been retrieved
            while record_count < max_total_records:
                # apply all the filters using the "and" syntax: [[filter1, filter2]]
                try:
                    csw.getrecords2(constraints=[fes_filters],
                                     esn='full',
                                         #outputschema="http://www.opengis.net/cat/csw/2.0.2",
                                         maxrecords=max_query_records,
                                         startposition=startposition)
                    logger.debug('CSW request:\n{}'.format(csw.request))
                    logger.debug('CSW response:\n {}'.format(csw.response))
                except Exception as e:
                    logger.warning('CSW request failed: {}'.format(e))
                    logger.debug('Bad CSW request:\n{}'.format(csw.request))
                    break
                    
    
                query_record_count = len(csw.records)
    
                for uuid in [uuid for uuid in csw.records.keys() if uuid not in uuid_list]:
                    record = csw.records[uuid]
                    title = record.title
    
                    #===========================================================
                    # # Ignore datasets with no distributions
                    # if not record.uris:
                    #     #logger.warning('No distribution(s) found for "{}"'.format(title)
                    #     continue
                    #===========================================================
    
                    uuid_list.append(uuid) # Remember UUID to avoid returning duplicates
                    
    #                logger.debug('bbox = {}'.format(record.bbox.__dict__)
                    #pprint(record.__dict__)
                    record_dict = {'csw': csw.url,
                                   'uuid': uuid,
                                   'title': title,
                                   'publisher': record.publisher,
                                   'author': record.creator,
                                   'abstract': record.abstract,
                                   'type': record.type,
                                  }
    
                    if record.bbox:
                        record_dict['bbox'] = [record.bbox.minx, record.bbox.minx, record.bbox.maxx, record.bbox.maxy],
                        record_dict['bbox_crs'] = record.bbox.crs or 'EPSG:4326'
    
                    # Deal with weird OWSLib behaviour where a single dict containing 'None' string values is returned
                    # when no distributions exist
                    if (len(record.uris) > 1 or 
                        ((len(record.uris) == 1) and record.uris[0]['url'] and (record.uris[0]['url'] != 'None'))
                        ):
                        distribution_info_list = copy.deepcopy(record.uris)
                    else:
                        # Create single dummy distribution for no protocol
                        distribution_info_list = [{'description': '',
                                                   'name': '',
                                                   'protocol': '',
                                                   'url': ''
                                                   }
                                                  ]
    
                    # Add layer information for web services
                    if get_layers:
                        for distribution_info in [distribution_info
                                                  for distribution_info in distribution_info_list
                                                  if distribution_info['protocol'] == 'OGC:WMS'
                                                  ]:
                            try:
                                wms = WebMapService(distribution_info['url'], version='1.1.1')
                                distribution_info['layers'] = wms.contents.keys()
                            except:
                                distribution_info['layers'] = ''
                                
        
                        for distribution_info in [distribution_info
                                                  for distribution_info in distribution_info_list
                                                  if distribution_info['protocol'] == 'OGC:WCS'
                                                  ]:
                            try:
                                wcs = WebCoverageService(distribution_info['url'], version='1.0.0')
                                distribution_info['layers'] = wcs.contents.keys()
                            except:
                                distribution_info['layers'] = ''
    
                    record_dict['distributions'] = distribution_info_list
                    record_dict['keywords'] = record.subjects
    
                    record_count += 1
                    yield record_dict
                    logger.debug('{} distribution(s) found for "{}"'.format(len(distribution_info_list), title))
    
                    if record_count >= max_total_records:  # Don't go around again for another query - maximum retrieved
                        raise Exception('Maximum number of records retrieved ({})'.format(max_total_records))
        
                if query_record_count < max_query_records:  # Don't go around again for another query - should be the end
                    break
    
                # Increment start position and repeat query
                startposition += max_query_records
    
        logger.debug('{} records found.'.format(record_count))


    def query_csw(self,
                  keyword_list=None,
                  bounding_box=None,
                  bounding_box_crs=None,
                  anytext_list=None,
                  titleword_list=None,
                  start_datetime=None,
                  stop_datetime=None,
                  record_type_list=None,
                  max_total_records=None,
                  get_layers=None
                  ):
        '''
        Function to query CSW using AND combination of provided search parameters and return generator object
            yielding nested dicts containing information about each record including distributions
        @param keyword_list: List of strings or comma-separated string containing keyword search terms
        @param bounding_box: Bounding box to search as a list of ordinates [bbox.minx, bbox.minx, bbox.maxx, bbox.maxy]
        @param bounding_box_crs: Coordinate reference system for bounding box. Defaults to value of self.settings['DEFAULT_CRS']
        @param anytext_list: List of strings or comma-separated string containing any text search terms
        @param titleword: List of strings or comma-separated string containing title search terms
        @param start_datetime: Datetime object defining start of temporal search period
        @param stop_datetime: Datetime object defining end of temporal search period
        @param record_type_list: List of strings or comma-separated string containing record type(s) to return
        @param max_total_records: Maximum total number of records to return. Defaults to value of self.settings['DEFAULT_MAXTOTALRECORDS']
        @param get_layers: Boolean flag indicating whether to get WMS/WCS layer names. Defaults to value of self.settings['DEFAULT_GET_LAYERS']
        
        @return: generator object yielding nested dicts containing information about each record including distributions
        '''
        bounding_box_crs = bounding_box_crs or self.settings['DEFAULT_CRS']
        get_layers = self.settings['DEFAULT_GET_LAYERS'] if get_layers is None else get_layers

        # Convert strings to lists if required
        if type(keyword_list) == str:
            keyword_list = self.list_from_comma_separated_string(keyword_list)

        if type(anytext_list) == str:
            anytext_list = self.list_from_comma_separated_string(anytext_list)

        if type(titleword_list) == str:
            titleword_list = self.list_from_comma_separated_string(titleword_list)

        record_type_list = record_type_list or self.settings['DEFAULT_RECORD_TYPES']
        if type(record_type_list) == str:
            record_type_list = self.list_from_comma_separated_string(record_type_list)

        # Build filter list
        fes_filter_list = []
        
        # Check for unchanged, upper-case, lower-case and capitalised keywords 
        # with single-character wildcards substituted for whitespace characters
        # GeoNetwork keyword search is always case sensitive
        if keyword_list:
            fes_filter_list += self.any_case_match_filters('Subject', keyword_list)
          
        if anytext_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='anyText', literal=phrase, matchCase=False) for phrase in anytext_list]
        
        if start_datetime or stop_datetime:
            fes_filter_list += self.get_date_filter(start_datetime, stop_datetime)
        
        if titleword_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='title', literal=titleword, matchCase=False) for titleword in titleword_list]
        
        # Check for unchanged, upper-case, lower-case and capitalised keywords 
        # with single-character wildcards substituted for whitespace characters
        # GeoNetwork type search is always case sensitive
        if record_type_list:
            fes_filter_list += self.any_case_match_filters('type', record_type_list)
        
        if bounding_box:
            # N.B: Bounding box ordinate ordering must match CRS. Default CRS84 supports lon-lat ordering, not lat-lon
            # See https://gis.stackexchange.com/questions/124050/how-do-i-specify-the-lon-lat-ordering-in-csw-bounding-box-request
            fes_filter_list += [fes.BBox(bounding_box, crs=bounding_box_crs)]

        assert fes_filter_list, 'No search criteria defined'

        # Use single filter if no "and" required
        if len(fes_filter_list) == 1:
            fes_filter_list = fes_filter_list[0]

        # Return generator object
        return self.get_csw_records(fes_filter_list, 
                                    max_total_records=max_total_records,
                                    get_layers=get_layers
                                    )

    def flatten_distribution_dict(self, record_dict, distribution_dict):
        '''
        Helper function to create a flattened dict from a record dict and a single distribution dict as returned by query_csw
        @param record_dict:
        @param distribution_dict
        
        @return dict: flattened dict containing information about specified record and distribution
        '''
        dataset_distribution_dict = copy.copy(record_dict) # Create shallow copy of record dict

        # Delete list of all distributions from copy of record dict
        del dataset_distribution_dict['distributions']

        # Convert lists to strings
        dataset_distribution_dict['keywords'] = ', '.join(sorted(dataset_distribution_dict['keywords']))
        
        bbox = dataset_distribution_dict.get('bbox')
        if bbox:
            dataset_distribution_dict['bbox'] = ', '.join(bbox[0]) #TODO: Cater for multiple bounding boxes

        # Merge distribution info into copy of record dict
        if distribution_dict:
            dataset_distribution_dict.update(distribution_dict)
            # Remove any leading " file://" from URL to give plain filename
            dataset_distribution_dict['url'] = re.sub('^file://', '', dataset_distribution_dict['url'])
            
            if 'layers' in dataset_distribution_dict.keys():
                dataset_distribution_dict['layers'] = ', '.join(dataset_distribution_dict['layers'])
        
        return dataset_distribution_dict
    
        
    def partial_string_match(self, superstring, partial_string_list):
        '''
        Helper function to perform partial string match for strings in a list
        @param superstring: Whole string against which to search for partial matches
        @param partial_string_list: List of partial strings for which to search superstring
        
        @return bool: True if any partial string found in superstring
        '''
        # Treat empty list as wildcard match
        if partial_string_list is None:
            return True
        
        for search_string in partial_string_list:
            if (not (search_string or superstring) # Both empty strings
                or (search_string and search_string in superstring) # Substring
                ):
                return True
        return False
            
            
    def get_distributions(self, distribution_protocols, dataset_dict_generator):
        '''
        Generator to yield flattened dicts containing information for all distributions matching 
        specified distribution_protocols
        @param distribution_protocols: comma-separated string or list containing distribution_protocols to 
        match (case insensitive partial string match). None or empty list treated as wildcard.
        @param dataset_dict_generator: Generator yeilding dict objects containing information about each record including distributions

        '''
        # Ensure protocol_list contains lower case strings for case insensitivity
        if type(distribution_protocols) == str:
            search_protocol_list = [distribution_protocol.lower() 
                                    for distribution_protocol in 
                                    self.list_from_comma_separated_string(distribution_protocols)]
        elif type(distribution_protocols) == list:
            search_protocol_list = [distribution_protocol.lower() 
                                    for distribution_protocol in distribution_protocols]
        else: 
            search_protocol_list = distribution_protocols
        
        for record_dict in dataset_dict_generator:
            if record_dict['distributions']:
                for distribution_dict in record_dict['distributions']:
                    # Replace None with empty string for all values
                    for key, value in iter(distribution_dict.items()):
                        if value is None:
                            distribution_dict[key] = ''
                                    
                    # If protocol match is found (case insensitive partial match) or 
                    if self.partial_string_match(distribution_dict['protocol'].lower(), search_protocol_list): 
                        yield self.flatten_distribution_dict(record_dict, distribution_dict)
            else: # This shouldn't be needed - "fake" distribution now created for no distributions
                yield self.flatten_distribution_dict(record_dict, None)
                

    def get_netcdf_urls(self, dataset_dict_generator):
        '''
        Generator to yield flattened dicts containing information for any netCDF distributions (file or OPeNDAP URL, file by preference)
        @param dataset_dict_generator: Generator yeilding dict objects containing information about each record including distributions
        '''
        for record_dict in dataset_dict_generator:
            distribution_dict = None
            for file_distribution in record_dict['distributions']:
                if 'file' in file_distribution['protocol'].lower():
                    match = re.match('(^file://)*(.*\.nc)$', file_distribution['url'])
                    try:
                        file_distribution['url'] = match.group(2) # Ignore any leading "file://"
                        if os.path.isfile(file_distribution['url']) and netCDF4.Dataset(file_distribution['url']): # Test for valid netCDF file
                            #logger.debug('file found')
                            distribution_dict = file_distribution
                            break
                    except:
                        logger.warning('Unable to open netCDF file {}'.format(file_distribution['url']))
                    
            if not distribution_dict:
                # Check for valid OPeNDAP endpoint if no valid file found
                for opendap_distribution in record_dict['distributions']:
                    if 'opendap' in opendap_distribution['protocol'].lower():
                        match = re.match('(.*\.nc)(\.html)*$', opendap_distribution['url'])
                        try:
                            opendap_distribution['url'] = match.group(1) # Ignore any trailing ".html"
                            if netCDF4.Dataset(opendap_distribution['url']): # Test for valid OPeNDAP endpoint
                                distribution_dict = opendap_distribution
                                break
                        except:
                            logger.warning('Unable to open OPeNDAP URL {}'.format(opendap_distribution['url']))
            
            if distribution_dict:
                yield self.flatten_distribution_dict(record_dict, distribution_dict)
                continue

            
    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, debug_value):
        if self._debug != debug_value or self._debug is None:
            self._debug = debug_value
            
            if self._debug:
                logger.setLevel(logging.DEBUG)
                logging.getLogger(self.__module__).setLevel(logging.DEBUG)
            else:
                logger.setLevel(logging.INFO)
                logging.getLogger(self.__module__).setLevel(logging.INFO)
                
        logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))
        logging.getLogger(self.__module__).debug('Logger {} set to level {}'.format(self.__module__, logger.level))

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
