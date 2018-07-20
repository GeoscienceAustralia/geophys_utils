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

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # Initial logging level for this module

DEBUG = False

DATABASE_ENGINE = 'SQLite'
#DATABASE_ENGINE = 'Postgres'

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
                    nc_file_template=None,
                    opendap_path_map=None
                    ):
        '''
        Function to populate DB with metadata from netCDF files
        '''
        nc_file_template = nc_file_template or '*.nc'
        
        for nc_path in self.find_files(nc_root_dir, file_template=nc_file_template):
            logger.info('Reading attributes from {}'.format(nc_path))
            nc_dataset = netCDF4.Dataset(nc_path, 'r')
            
            nc_attribute = dict(nc_dataset.__dict__)
    
            #nc_attribute['nc_path'] = nc_path
            nc_attribute['nc_path'] = nc_path.replace('D:\\Temp\\gravity point data\\', '/g/data2/uc0/rr2_dev/axi547/ground_gravity/point_datasets/') #TODO: Remove this temporary hack
            
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
            
            if opendap_path_map:
                distribution_list.append(Distribution(url=nc_attribute['nc_path'].replace(*opendap_path_map),
                                                      protocol='opendap'
                                                      )
                                         )
                
            dataset = Dataset(dataset_title=nc_attribute['title'],
                              ga_survey_id=nc_attribute['survey_id'],
                              longitude_min=np.asscalar(nc_attribute['geospatial_lon_min']),
                              longitude_max=np.asscalar(nc_attribute['geospatial_lon_max']),
                              latitude_min=np.asscalar(nc_attribute['geospatial_lat_min']),
                              latitude_max=np.asscalar(nc_attribute['geospatial_lat_max']),
                              convex_hull_polygon=nc_attribute.get('geospatial_bounds'), 
                              keyword_list=[keyword.strip() for keyword in nc_attribute['keywords'].split(',')],
                              distribution_list=distribution_list,
                              metadata_uuid=nc_attribute.get('uuid') # Could be None
                              )

            #logger.debug('dataset: {}'.format(dataset.__dict__))
            self.dataset_metadata_cache.add_dataset(dataset)
        

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
                       opendap_path_map=('/g/data2/uc0/rr2_dev/',
                                         'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/'
                                         )
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
