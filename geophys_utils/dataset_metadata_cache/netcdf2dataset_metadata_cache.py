'''
Created on 20 Jul. 2018

@author: Alex
'''
import sys
import os
from glob import glob
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache, Dataset, Distribution
import logging
import netCDF4
import numpy as np
from datetime import datetime, date

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

DEBUG = True

DATABASE_ENGINE = 'SQLite'
#DATABASE_ENGINE = 'Postgres'

# Set this to change file path - THIS IS A HACK!!!!
#FILE_PATH_MAPS = None
FILE_PATH_MAPS = {'D:\\Temp\\gravity point_datasets\\': '/g/data2/uc0/rr2_dev/axi547/ground_gravity/point_datasets/',
                  'D:\\Temp\\AEM conductivity datasets\\': '/g/data2/uc0/rr2_dev/axi547/aem/'
                  }

# Set this to derive OPeNDAP endpoint URL from file path - THIS IS A HACK!!!!
#OPENDAP_PATH_MAPS = None
OPENDAP_PATH_MAPS = {'/g/data2/uc0/rr2_dev/': 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/',
                     '/g/data1/rr2/': 'http://dapds00.nci.org.au/thredds/dodsC/rr2/'
                     }

class NetCDF2DatasetMetadataCache(object):
    '''
    classdocs
    '''


    def __init__(self, debug=True):
        '''
        Constructor
        '''
        self.dataset_metadata_cache = get_dataset_metadata_cache(db_engine=DATABASE_ENGINE, debug=debug)
        
        
    def find_files(self, root_dir, file_template, extension_filter='.nc'):
        '''
        Function to simulate the result of a filtered Linux find command
        Uses glob with user-friendly file system wildcards instead of regular expressions for template matching
        @param root_dir: Top level directory to be searched
        @param file_template: glob-style filename template similar to -name argument of linux find command
        @param extension_filter: single file extension (including ".") on which to filter files
        @return file_path_list: List of file paths
        '''
        #===========================================================================
        # file_path_list = sorted([filename for filename in subprocess.check_output(
        #     ['find', args.netcdf_dir, '-name', args.file_template]).split('\n') if re.search('\.nc$', filename)])
        #===========================================================================
        root_dir = os.path.abspath(root_dir)
        file_path_list = glob(os.path.join(root_dir, file_template))
        for topdir, subdirs, _files in os.walk(root_dir, topdown=True):
            for subdir in subdirs:
                file_path_list += [file_path 
                                   for file_path in glob(os.path.join(topdir, subdir, file_template))
                                   if os.path.isfile(file_path)
                                   and os.path.splitext(file_path)[1] == extension_filter
                                   ]
        file_path_list = sorted(file_path_list)    
        return file_path_list

    def populate_db(self,
                    nc_root_dir,
                    nc_file_template=None
                    ):
        '''
        Function to populate DB with metadata from netCDF files
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
        
            
        nc_file_template = nc_file_template or '*.nc'
        
        for nc_path in self.find_files(nc_root_dir, file_template=nc_file_template):
            logger.info('Reading attributes from {}'.format(nc_path))
            nc_dataset = netCDF4.Dataset(nc_path, 'r')
            
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
            distribution_list = [Distribution(url='file://'+nc_attribute['nc_path'],
                                              protocol='file'
                                              )
                                 ]
            
            if OPENDAP_PATH_MAPS:
                # Replace directory with URL prefix
                opendap_url=nc_attribute['nc_path']
                for opendap_path_map in OPENDAP_PATH_MAPS.items():
                    opendap_url = opendap_url.replace(*opendap_path_map)
                    
                distribution_list.append(Distribution(url=opendap_url,
                                                      protocol='opendap'
                                                      )
                                         )
                
            try:
                point_count = nc_dataset.dimensions['point'].size
            except:
                point_count = None
            
            try:
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
    assert len(sys.argv) >= 2 and len(sys.argv) <= 3, 'Usage: {} <root_dir> [<file_template>]'.format(sys.argv[0])
    nc_root_dir = sys.argv[1]
    if len(sys.argv) == 3:
        nc_file_template = sys.argv[2]
    else:
        nc_file_template = '*.nc'
    

    nc2dmc = NetCDF2DatasetMetadataCache(debug=DEBUG)
    nc2dmc.populate_db(nc_root_dir,
                       nc_file_template=nc_file_template,
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
