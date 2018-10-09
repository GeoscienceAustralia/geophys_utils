'''
Created on 9 Oct. 2018

@author: Alex Ip

Utility to read metadata from CSW query and/or netCDF file into dataset metadata cache
'''
import sys
import re
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache, Dataset, Distribution
import logging
import netCDF4
import numpy as np
from datetime import datetime
from geophys_utils import CSWUtils

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

DEBUG = True

DATABASE_ENGINE = 'SQLite'
#DATABASE_ENGINE = 'Postgres'

DEFAULT_CSW_URL = 'https://ecat.ga.gov.au/geonetwork/srv/eng/csw'

# Set this to derive file path from OPeNDAP endpoint URL - THIS IS A HACK!!!!
#OPENDAP_PATH_MAPS = None
FILE_PATH_MAPS = {'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/': '/g/data2/uc0/rr2_dev/',
                     'http://dapds00.nci.org.au/thredds/dodsC/rr2/': '/g/data1/rr2/'
                     }

class CSW2DatasetMetadataCache(object):
    '''
    classdocs
    '''


    def __init__(self, debug=True):
        '''
        Constructor
        '''
        self.dataset_metadata_cache = get_dataset_metadata_cache(db_engine=DATABASE_ENGINE, debug=debug)
        
        
    def populate_db(self,
                    keyword_list=None,
                    anytext=None,
                    titleword_list=None,
                    bounding_box=None,
                    start_datetime=None,
                    stop_datetime=None,
                    record_type_list=None,
                    csw_url=None
                    ):
        '''
        Function to populate DB with metadata from CSW query and/or netCDF files via OPeNDAP
        '''
        
        def datetimestring2date(datetime_string):
            '''
            Function to return Python date object from a datetime string
            '''
            result = None
            
            #logger.debug('datetime_string: {}'.format(datetime_string))
            
            if datetime_string:
                for datetime_format in ['%Y-%m-%d', 
                                        '%Y-%m-%d %H:%M:%S'
                                        ]:
                    try:
                        result = datetime.strptime(datetime_string, datetime_format).date()
                        break
                    except ValueError:
                        pass
            
            return result
        
        csw_url = csw_url or DEFAULT_CSW_URL
        csw_utils = CSWUtils(csw_url_list=[csw_url], 
                             timeout=None,
                             debug=DEBUG,
                             settings_path=None
                             )
        
        record_generator = csw_utils.query_csw(keyword_list=keyword_list,
                                     anytext_list=anytext,
                                     titleword_list=titleword_list,
                                     bounding_box=bounding_box,
                                     start_datetime=start_datetime,
                                     stop_datetime=stop_datetime,
                                     record_type_list=record_type_list,
                                     max_total_records=None,
                                     get_layers=None
                                     )
        
        distribution_count = 0
        for distribution_dict in csw_utils.get_distributions(['opendap'], record_generator):
            distribution_count += 1
            
            opendap_url = re.sub('\.html$', '', distribution_dict['url']) # Strip trailing ".html"
            
            # Derive NCI file path from OPeNDAP URL
            nc_path = opendap_url
            if FILE_PATH_MAPS:
                for file_path_map in FILE_PATH_MAPS.items():
                    nc_path = nc_path.replace(*file_path_map)
        
            logger.info('Reading attributes from {}'.format(opendap_url))
            nc_dataset = netCDF4.Dataset(opendap_url, 'r')
             
            nc_attribute = dict(nc_dataset.__dict__)
     
            if FILE_PATH_MAPS:
                for file_path_map in FILE_PATH_MAPS.items():
                    nc_path = nc_path.replace(*file_path_map)
                 
            nc_attribute['nc_path'] = nc_path
            
#===============================================================================
# // global attributes:
#                 :_NCProperties = "version=1|netcdflibversion=4.5.0|hdf5libversion=1.10.0" ;
#                 :title = "Batten Fault Zone Gravity Survey" ;
#                 :survey_id = "201780" ;
#                 :Conventions = "CF-1.6,ACDD-1.3" ;
#                 :keywords = "points, gravity, ground digital data, geophysical survey, survey 201780, AUS, NT, Earth sciences, geophysics, geoscientificInformation" ;
#                 :geospatial_lon_min = 135.3741f ;
#                 :geospatial_lon_max = 137.1358f ;
#                 :geospatial_lon_units = "degrees_east" ;
#                 :geospatial_long_resolution = "point" ;
#                 :geospatial_lat_min = -18.23301f ;
#                 :geospatial_lat_max = -15.33811f ;
#                 :geospatial_lat_units = "degrees_north" ;
#                 :geospatial_lat_resolution = "point" ;
#                 :history = "Pulled from point gravity database at Geoscience Australia" ;
#                 :summary = "This gravity survey, 201780, Batten Fault Zone Gravity Survey located in NT measures the slight variations in the earth\'s gravity based on the underlying structure or geology" ;
#                 :location_accuracy_min = 1.f ;
#                 :location_accuracy_max = 1.f ;
#                 :time_coverage_start = "2017-09-28 00:00:00" ;
#                 :time_coverage_end = "2017-11-26 00:00:00" ;
#                 :time_coverage_duration = "59 days, 0:00:00" ;
#                 :date_created = "2018-07-14T02:34:58.573387" ;
#                 :institution = "Geoscience Australia" ;
#                 :source = "ground observation" ;
#                 :cdm_data_type = "Point" ;
#                 :geospatial_bounds = "POLYGON((135.5491 -18.2330, 135.3972 -18.2330, 135.3786 -18.2151, 135.3766 -18.1796, 135.3741 -17.1306, 135.4804 -15.4322, 135.5047 -15.4098, 135.5364 -15.3999, 135.6383 -15.3767, 135.7226 -15.3602, 135.8378 -15.3424, 136.0260 -15.3381, 136.0985 -15.3559, 137.1116 -15.9800, 137.1358 -15.9956, 137.1015 -18.2200, 137.0246 -18.2236, 136.6086 -18.2277, 136.4190 -18.2295, 135.5491 -18.2330))" ;
#                 :_Format = "netCDF-4" ;
# }
#===============================================================================
            distribution_list = [Distribution(url=opendap_url,
                                              protocol='opendap'
                                              ),
                                 Distribution(url='file://'+nc_attribute['nc_path'],
                                              protocol='file'
                                              )            
                                 
                                 ]                
            try:
                point_count = nc_dataset.dimensions['point'].size
            except:
                point_count = None
            
            try:
                #TODO: Read stuff from CSW result as well as netCDF dataset
                dataset = Dataset(dataset_title=nc_attribute['title'],
                                  ga_survey_id=nc_attribute.get('survey_id'),
                                  longitude_min=np.asscalar(nc_attribute['geospatial_lon_min']),
                                  longitude_max=np.asscalar(nc_attribute['geospatial_lon_max']),
                                  latitude_min=np.asscalar(nc_attribute['geospatial_lat_min']),
                                  latitude_max=np.asscalar(nc_attribute['geospatial_lat_max']),
                                  convex_hull_polygon=nc_attribute.get('geospatial_bounds'), 
                                  keyword_list=[keyword.strip() for keyword in nc_attribute['keywords'].split(',')],
                                  distribution_list=distribution_list,
                                  point_count=point_count,
                                  metadata_uuid=nc_attribute.get('uuid'), # Could be None
                                  start_date=datetimestring2date(nc_attribute.get('time_coverage_start')),
                                  end_date=datetimestring2date(nc_attribute.get('time_coverage_end'))
                                  )

                #logger.debug('dataset: {}'.format(dataset.__dict__))
                self.dataset_metadata_cache.add_dataset(dataset)
            except Exception as e:
                logger.warning('Unable to process dataset {}: {}'.format(nc_path, e))
        

def main():
    assert len(sys.argv) >= 2 and len(sys.argv) <= 3, 'Usage: {} <keyword_list>'.format(sys.argv[0])
    keyword_list = sys.argv[1] # Should be comma-separated keyword list
    #TODO: Allow for more command line options for CSW search (copy stuff from _csw_utils.py?)

    csw2dmc = CSW2DatasetMetadataCache(debug=DEBUG)
    csw2dmc.populate_db(keyword_list=keyword_list
                       )
    
        
if __name__ == '__main__':
    # Setup logging handlers if required
    if not logger.handlers:
        # Set handler for root logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        #console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for logger {}'.format(logger.name))

    main()
