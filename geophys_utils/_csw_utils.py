'''
Created on 23Feb.,2017

@author: u76345
'''
import re
import copy
from pprint import pprint
from owslib import fes
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
    DEFAULT_MAXTOTALRECORDS = 500 # Maximum total number of records to retrieve

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

    def get_csw_info(self, fes_filters, maxrecords=None):
        '''
        Function to find all distributions for all records returned
        Returns a nested dict keyed by UUID
        @param fes_filters: List of fes filters to apply to CSW query
        @param maxrecords: Maximum number of records to return per CSW query. Defaults to value of CSWUtils.DEFAULT_MAXRECORDS
        
        @return: Nested dict object containing information about each record including distributions
        '''
        dataset_dict = {} # Dataset details keyed by title
    
        maxrecords = maxrecords or CSWUtils.DEFAULT_MAXRECORDS 
        startposition = 1 # N.B: This is 1-based, not 0-based
        
        while True: # Keep querying until all results have been retrieved
            # apply all the filters using the "and" syntax: [[filter1, filter2]]
            self.csw.getrecords2(constraints=[fes_filters], 
                                 esn='full', 
                                 maxrecords=maxrecords, 
                                 startposition=startposition)
            
    #        print 'csw.request = %s' % csw.request
    #        print 'csw.response = %s' % csw.response
            
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

            if len(dataset_dict) >= CSWUtils.DEFAULT_MAXTOTALRECORDS:  # Don't go around again for another query - maximum retrieved
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
                  stop_datetime=None
                  ):
        '''
        Function to query CSW using AND combination of provided search parameters
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
        if keyword_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='Subject', literal=keyword, matchCase=False) for keyword in keyword_list]
        if anytext_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='anyText', literal=phrase, matchCase=False) for phrase in anytext_list]
        if start_datetime or stop_datetime:
            fes_filter_list += self.get_date_filter(start_datetime, stop_datetime)
        if titleword_list:
            fes_filter_list += [fes.PropertyIsLike(propertyname='title', literal=titleword, matchCase=False) for titleword in titleword_list]
        if bounding_box:
            fes_filter_list += [fes.BBox(bounding_box, crs=bounding_box_crs)]


        if len(fes_filter_list) == 1:
            fes_filter_list = fes_filter_list[0]
        print fes_filter_list

        return self.get_csw_info(fes_filter_list)
            
def main():
    '''
    Quick and dirty main function for on-the-fly testing
    '''
    #keywords = "Geophysical National Coverage, NCI, geoscientific%Information, grid" # National Coverages
    keywords = "TMI, magnetics, NCI, AU, Magnetism and Palaeomagnetism, Airborne Digital Data, Geophysical Survey, grid" # Magnetic survey grids
    #titlewords = "onshore, gravity, grid, Australia, 2016"
    anytext = "Magnetism, Palaeomagnetism"
    bounds = [148.996, -35.48, 149.399, -35.124]

    cswu = CSWUtils()
    
    result_dict = cswu.query_csw(keyword_list=keywords, bounding_box=bounds, anytext_list=anytext)
    #result_dict = cswu.query_csw(titleword_list=titlewords)


    pprint(result_dict)
    
    print '%d results found.' % len(result_dict)

if __name__ == '__main__':
    main()