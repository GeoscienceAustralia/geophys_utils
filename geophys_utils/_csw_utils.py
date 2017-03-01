'''
Created on 23Feb.,2017

@author: u76345
'''
import re
import copy
from pprint import pprint
from datetime import datetime, timedelta
from owslib import fes
import argparse
from owslib.csw import CatalogueServiceWeb
from owslib.wms import WebMapService
from owslib.wcs import WebCoverageService

class CSWUtils(object):
    '''
    CSW query utilities
    '''
    DEFAULT_CSW_URL = 'http://ecat.ga.gov.au/geonetwork/srv/eng/csw' # GA's externally-facing eCat
    DEFAULT_TIMEOUT = 30 # Timeout in seconds
    DEFAULT_CRS = 'EPSG:4326' # Unprojected WGS84
    DEFAULT_MAXRECORDS = 100 # Retrieve only this many datasets per CSW query
    DEFAULT_MAXTOTALRECORDS = 1000 # Maximum total number of records to retrieve

    def __init__(self, csw_url=None, timeout=None):
        '''
        Constructor for CSWUtils class
        @param csw_url: URL for CSW service. Defaults to value of CSWUtils.DEFAULT_CSW_URL
        @param timeout: Timeout in seconds. Defaults to value of CSWUtils.DEFAULT_TIMEOUT
        '''
        csw_url = csw_url or CSWUtils.DEFAULT_CSW_URL
        timeout = timeout or CSWUtils.DEFAULT_TIMEOUT

        self.csw = CatalogueServiceWeb(csw_url, timeout=timeout)

    def list_from_comma_separated_string(self, comma_separated_string):
        '''
        Helper function to return list of strings from a comma-separated string
        Substitute single-character wildcard for whitespace characters
        @param comma_separated_string: comma-separated string
        @return: list of strings
        '''
        return [re.sub('(\s)', '_', keyword.strip()) for keyword in comma_separated_string.split(',')]

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
            start_date_string = start_datetime.isoformat()
        else:
            start_date_string = '1900-01-01T00:00:00'

        if start_datetime:
            stop_date_string = start_datetime.isoformat()
        else:
            stop_date_string = '2100-01-01T23:59:59'

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
        Resolves issue where "matchCase=False" is ignored
        @param propertyname: String denoting CSW property to search
        @param word_list: List of strings for (partially) case-insensitive "and" search
        
        @return: List of FES filters expressing "and" query
        '''
        return [fes.Or([fes.PropertyIsLike(propertyname=propertyname, literal=word_variant, matchCase=False)
                        for word_variant in set([word, word.upper(), word.lower(), word.title()])
                        ]
                       ) for word in word_list
                ]

    def get_csw_info(self, fes_filters, maxrecords=None, max_total_records=None):
        '''
        Function to find all distributions for all records returned
        Returns a nested dict keyed by UUID
        @param fes_filters: List of fes filters to apply to CSW query
        @param maxrecords: Maximum number of records to return per CSW query. Defaults to value of CSWUtils.DEFAULT_MAXRECORDS
        @param max_total_records: Maximum total number of records to return. Defaults to value of CSWUtils.DEFAULT_MAXTOTALRECORDS

        @return: Nested dict object containing information about each record including distributions
        '''
        dataset_dict = {} # Dataset details keyed by title

        maxrecords = maxrecords or CSWUtils.DEFAULT_MAXRECORDS
        max_total_records = max_total_records or CSWUtils.DEFAULT_MAXTOTALRECORDS

        startposition = 1 # N.B: This is 1-based, not 0-based

        # Keep querying until all results have been retrieved
        while len(dataset_dict) < max_total_records:
            # apply all the filters using the "and" syntax: [[filter1, filter2]]
            self.csw.getrecords2(constraints=[fes_filters],
                                 esn='full',
                                 maxrecords=maxrecords,
                                 startposition=startposition)

#            print 'csw.request = %s' % self.csw.request
#            print 'self.csw.response = %s' % self.csw.response

            record_count = len(self.csw.records)

            for uuid in self.csw.records.keys():
                record = self.csw.records[uuid]
                title = record.title

                # Ignore datasets with no distributions
                if not record.uris:
                    #print 'No distribution(s) found for "%s"' % title
                    continue

#                print 'bbox = %s' % record.bbox.__dict__

                record_dict = {'uuid': uuid,
                               'title': title,
                               'publisher': record.publisher,
                               'author': record.creator,
                               'abstract': record.abstract,
                              }

                if record.bbox:
                    record_dict['bbox'] = [record.bbox.minx, record.bbox.minx, record.bbox.maxx, record.bbox.maxy],
                    record_dict['bbox_crs'] = record.bbox.crs or 'EPSG:4326'

                distribution_info_list = copy.deepcopy(record.uris)

                # Add layer information for web services
                for distribution_info in [distribution_info
                                          for distribution_info in distribution_info_list
                                          if distribution_info['protocol'] == 'OGC:WMS'
                                          ]:
                    wms = WebMapService(distribution_info['url'], version='1.1.1')
                    distribution_info['layers'] = wms.contents.keys()

                for distribution_info in [distribution_info
                                          for distribution_info in distribution_info_list
                                          if distribution_info['protocol'] == 'OGC:WCS'
                                          ]:
                    wcs = WebCoverageService(distribution_info['url'], version='1.0.0')
                    distribution_info['layers'] = wcs.contents.keys()

                record_dict['distributions'] = distribution_info_list
                record_dict['keywords'] = record.subjects

                dataset_dict[uuid] = record_dict
                #print '%d distribution(s) found for "%s"' % (len(info_list), title)

                if len(dataset_dict) >= max_total_records:  # Don't go around again for another query - maximum retrieved
                    break
    
            if record_count < maxrecords:  # Don't go around again for another query - should be the end
                break

            startposition += maxrecords

        #assert distribution_dict, 'No URLs found'
        #print '%d records found.' % len(dataset_dict)
        return dataset_dict

    def query_csw(self,
                  keyword_list=None,
                  bounding_box=None,
                  bounding_box_crs=None,
                  anytext_list=None,
                  titleword_list=None,
                  start_datetime=None,
                  stop_datetime=None,
                  max_total_records=None
                  ):
        '''
        Function to query CSW using AND combination of provided search parameters
        @param keyword_list: List of strings or comma-separated string containing keyword search terms
        @param bounding_box: Bounding box to search as a list of ordinates [bbox.minx, bbox.minx, bbox.maxx, bbox.maxy]
        @param bounding_box_crs: Coordinate reference system for bounding box. Defaults to value of CSWUtils.DEFAULT_CRS
        @param anytext_list: List of strings or comma-separated string containing any text search terms
        @param titleword: List of strings or comma-separated string containing title search terms
        @param start_datetime: Datetime object defining start of temporal search period
        @param stop_datetime: Datetime object defining end of temporal search period
        @param max_total_records: Maximum total number of records to return. Defaults to value of CSWUtils.DEFAULT_MAXTOTALRECORDS

        @return: Nested dict object containing information about each record including distributions
        '''
        bounding_box_crs = bounding_box_crs or CSWUtils.DEFAULT_CRS

        # Convert strings to lists if required
        if type(keyword_list) == str:
            keyword_list = self.list_from_comma_separated_string(keyword_list)

        if type(anytext_list) == str:
            anytext_list = self.list_from_comma_separated_string(anytext_list)

        if type(titleword_list) == str:
            titleword_list = self.list_from_comma_separated_string(titleword_list)

        # Build filter list
        fes_filter_list = []
        
        # Check for unchanged, upper-case, lower-case and capitalised keywords
        if keyword_list:
            fes_filter_list += self.any_case_match_filters('Subject', keyword_list)
          
        if anytext_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='anyText', literal=phrase, matchCase=False) for phrase in anytext_list]
        
        if start_datetime or stop_datetime:
            fes_filter_list += self.get_date_filter(start_datetime, stop_datetime)
        
        if titleword_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='title', literal=titleword, matchCase=False) for titleword in titleword_list]
        
        if bounding_box:
            fes_filter_list += [fes.BBox(bounding_box, crs=bounding_box_crs)]

        assert fes_filter_list, 'No search criteria defined'

        if len(fes_filter_list) == 1:
            fes_filter_list = fes_filter_list[0]

        return self.get_csw_info(fes_filter_list, 
                                 max_total_records=max_total_records)

    def find_distributions(self, distribution_protocol, dataset_dict):
        '''
        Function to return flattened list of dicts containing information for all
        distributions matching specified distribution_protocol
        @param distribution_protocol: distribution_protocol to match (case insensitive partial string match)
        @param dataset_dict: Nested dict object containing information about each record including distributions
        
        @return: flattened list of dicts containing information for all  distributions matching specified 
        distribution_protocol
        '''
        result_list = []
        for record_dict in dataset_dict.values():
            for distribution_dict in record_dict['distributions']:
                if distribution_protocol.upper() in distribution_dict['protocol'].upper(): # If protocol match is found
                    dataset_distribution_dict = copy.copy(record_dict) # Create shallow copy of record dict

                    # Delete list of all distributions from copy of record dict
                    del dataset_distribution_dict['distributions']

                    # Convert lists to strings
                    dataset_distribution_dict['keywords'] = ', '.join(dataset_distribution_dict['keywords'])
                    dataset_distribution_dict['bbox'] = ', '.join([str(ordinate) for ordinate in dataset_distribution_dict['bbox']])

                    # Merge distribution info into copy of record dict
                    dataset_distribution_dict.update(distribution_dict)
                    # Remove any leading " file://" from URL to give plain filename
                    dataset_distribution_dict['url'] = re.sub('^file://', '', dataset_distribution_dict['url'])

                    result_list.append(dataset_distribution_dict)

        return result_list


def date_string2datetime(date_string):
    '''
    Helper function to convert date string in one of several possible formats to a datetime object
    @param date_string: date string in one of several possible formats
    
    @return datetime object
    '''
    DATE_FORMAT_LIST = ['%Y%m%d', '%d/%m/%y', '%d/%m/%Y']
    if date_string:
        for format_string in DATE_FORMAT_LIST:
            try:
                datetime_result = datetime.strptime(date_string, format_string)
                break
            except ValueError:
                pass
            
        return datetime_result


def main():
    '''
    Main function to take command line parameters, perform CSW query and print required output
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-k", "--keywords", help="comma-separated list of keywords for search", type=str)
    parser.add_argument("-t", "--titlewords", help="comma-separated list of titlewords for search", type=str)
    parser.add_argument("-a", "--anytext", help="comma-separated list of text snippets for search", type=str)
    parser.add_argument("-b", "--bounds", help="comma-separated <minx>,<miny>,<maxx>,<maxy> ordinates of bounding box for search",
                        type=str)
    parser.add_argument("-c", "--crs", help='coordinate reference system for bounding box coordinates for search. Defaults to "EPSG:4326".',
                        type=str)
    parser.add_argument("-s", "--start_date", help="start date for search", type=str)
    parser.add_argument("-e", "--end_date", help="end date for search", type=str)
    # Parameters to define output
    parser.add_argument("-p", "--protocols", help='comma-separated list of distribution protocols for output. Defaults to "file".', type=str)
    parser.add_argument("-f", "--fields", help='comma-separated list of fields for output. Defaults to "protocol,url,title"', type=str)
    parser.add_argument("-d", "--delimiter", help='field delimiter for output. Defaults to "\t"', type=str)
    parser.add_argument("-u", "--url", help="CSW URL to query - Defaults to GA's external eCat", type=str)
    parser.add_argument("-m", "--max_results", help="Maximum number of records to return. Defaults to %d" % CSWUtils.DEFAULT_MAXTOTALRECORDS, type=int)
    parser.add_argument('-r', '--header_row', action='store_const', const=True, default=False,
                        help='display header row. Default is no header')
    
    args = parser.parse_args()

    # Convert string to list of floats
    if args.bounds:
        bounds = [float(ordinate) for ordinate in args.bounds.split(',')]
    else:
        bounds = None

    # Convert string to datetime
    start_date = date_string2datetime(args.start_date)
    #if args.start_date:        
        #print 'start_date = "%s"' % start_date.isoformat()

    # Convert string to datetime
    end_date = date_string2datetime(args.end_date)
    if args.end_date:
        end_date += timedelta(days=1) # Add one day to make end date inclusive        
        #print 'end_date = "%s"' % end_date.isoformat()
        
    # Default to listing file path
    protocol_list = ([protocol.strip() for protocol in args.protocols.split(',')] if args.protocols else None) or ['file']
    
    # Default to showing URL and title
    field_list = ([field.strip().lower() for field in args.fields.split(',')] if args.fields else None) or ['protocol', 'url', 'title']
    
    # Set default delimiter to tab character
    delimiter = args.delimiter or '\t'
    
    #create a CSW object and populate the parameters with argparse inputs - print results
    cswu = CSWUtils(args.url)
    result_dict = cswu.query_csw(keyword_list=args.keywords,
                                 anytext_list=args.anytext,
                                 titleword_list=args.titlewords,
                                 bounding_box=bounds,
                                 start_datetime=start_date,
                                 stop_datetime=end_date,
                                 max_total_records=args.max_results
                                 )
    #pprint(result_dict)
    #print '%d results found.' % len(result_dict)

    # Print header if required
    if args.header_row:
        print delimiter.join(['"' + field + '"' for field in field_list])
    
    # Print results
    for distribution_protocol in protocol_list:
        for distribution in cswu.find_distributions(distribution_protocol, result_dict):
            print delimiter.join(['"' + distribution[field] + '"' for field in field_list])


if __name__ == '__main__':
    main()
