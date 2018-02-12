'''
Created on 25 Jan. 2018

@author: Alex Ip

Quick-and-dirty point fetching demo
'''
import os
import sys
import csv
import netCDF4
import argparse
import logging
import re
from datetime import datetime
from pprint import pformat
import numpy as np

from geophys_utils import NetCDFLineUtils

# Setup logging handlers if required
logger = logging.getLogger(__name__) # Get __main__ logger
logger.setLevel(logging.INFO) # Initial logging level for this module
    
if not logger.handlers:
    # Set handler for root logger to standard output
    console_handler = logging.StreamHandler(sys.stdout)
    #console_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

class GeophysPointFetcher(object):
    '''
    GeophysPointFetcher class definition
    '''
    DEFAULT_METADATA_CSV_PATH = os.path.join(os.path.dirname(__file__), 'geophysics_line_nc_metadata.csv')
    OPENDAP_PATH_MAP = ('/g/data2/uc0', 'http://dapds00.nci.org.au/thredds/dodsC/uc0')

    def __init__(self, metadata_csv_path=None):
        '''
        Constructor
        '''
        metadata_csv_path = metadata_csv_path or GeophysPointFetcher.DEFAULT_METADATA_CSV_PATH
        
        #TODO: Replace this bodgy CSV-based code with a CSW catalogue based solution
        logger.debug('Importing metadata CSV file {}'.format(metadata_csv_path))
        with open(metadata_csv_path) as metadata_csv_file:
            self._metadata_keys = [key.lower() for key in csv.DictReader(metadata_csv_file).fieldnames]

            csv_reader = csv.reader(metadata_csv_file)
            
            self._metadata = [dict(zip(self._metadata_keys, [None if value == '' else value 
                                                             for value in row
                                                             ]
                                       )
                                   )
                              for row in csv_reader
                              ]    
            
        metadata_csv_file.close()
        
        
    def dataset_metadata_generator(self, bounding_box, keywords=[], metadata_filter_function=None):
        '''
        Generator returning all metadata dicts near given coordinate
        '''
        keyword_set = set(keywords)      
        
        #TODO: Replace this bodgy CSV-based code with a CSW catalogue based solution
        for metadata_dict in self._metadata:
            # Exclude any datasets not touching bounding box
            try:
                if not (float(metadata_dict['geospatial_lon_min']) <= bounding_box[0]
                    and float(metadata_dict['geospatial_lat_min']) <= bounding_box[1]
                    and float(metadata_dict['geospatial_lon_max']) >= bounding_box[2]
                    and float(metadata_dict['geospatial_lat_max']) >= bounding_box[3]
                    ):
                    continue
            except KeyError:
                logger.warning('Unable to determine bounding box for {}'.format(metadata_dict['file_path']))
                logger.debug('metadata_dict = {}'.format(pformat(metadata_dict)))
                continue
            
            # Exclude any datasets which do not contain all specified keywords
            if keyword_set and not keyword_set <= set([keyword.strip() for keyword in metadata_dict['keywords'].split(',')]):
                continue
        
            # Yield only dataset metadata satisfying filter function
            if not metadata_filter_function or metadata_filter_function(metadata_dict):
                yield metadata_dict


    def point_data_generator(self, 
                             bounding_box, 
                             keywords=[], 
                             metadata_filter_function=None, 
                             variable_names=None, 
                             flight_lines_only=True
                             ):
        '''
        Generator yielding point data for each dataset
        '''
        t0 = datetime.now()
        for metadata_dict in self.dataset_metadata_generator(bounding_box, 
                                                            keywords=keywords,
                                                            metadata_filter_function=metadata_filter_function
                                                            ):
            nc_path = metadata_dict['file_path']
            
            if not os.path.isfile(nc_path): 
                nc_path = nc_path.replace(GeophysPointFetcher.OPENDAP_PATH_MAP[0], 
                                          GeophysPointFetcher.OPENDAP_PATH_MAP[1]
                                          )
                
            try:
                logger.info('Opening {}'.format(nc_path))
                nc_dataset = netCDF4.Dataset(nc_path, 'r')
                netcdf_line_utils = NetCDFLineUtils(nc_dataset)
                
                if flight_lines_only:
                    print 'Excluding tie-lines'
                    line_numbers = nc_dataset.variables['line'][nc_dataset.variables['flag_linetype'][:] == 2]
                    line_mask = np.zeros(shape=(netcdf_line_utils.point_count,), dtype=bool)
                    for _line_number, single_line_mask in netcdf_line_utils.get_line_masks(line_numbers):
                        line_mask = np.logical_or(line_mask, single_line_mask)
                else:
                    line_mask = np.ones(shape=(netcdf_line_utils.point_count,).shape, dtype=bool)
                
                print 'Computing spatial subset mask'
                spatial_mask = netcdf_line_utils.get_spatial_mask(bounding_box)
                if not np.any(spatial_mask):
                    logger.warning('No points in bounding box {} for {}'.format(tuple(bounding_box), nc_path))
                    continue
                
                point_indices = np.where(np.logical_and(spatial_mask,
                                                        line_mask
                                                        )
                                         )[0]
            
                if not len(point_indices):
                    logger.warning('No points in bounding box {} for {}'.format(tuple(bounding_box), nc_path))
                    continue
                
                logger.info('{} points found in bounding box {} for {}'.format(len(point_indices), tuple(bounding_box), nc_path))

                dataset_dict = {}
                
                dataset_dict['metadata'] = dict(metadata_dict) # Copy metadata
                logger.debug('\tReading coordinates')
                dataset_dict['coordinates'] = netcdf_line_utils.xycoords[point_indices]
                dataset_dict['values'] = {}
                for variable_name in netcdf_line_utils.point_variables:
                    if not variable_names or variable_name in variable_names:
                        logger.debug('\tReading values for {}'.format(variable_name))
                        dataset_dict['values'][variable_name] = nc_dataset.variables[variable_name][:][point_indices]
                
                yield dataset_dict
            except Exception as e:
                logger.error('Failed to retrieve point data from {}: {}'.format(nc_path, e.message))
                
        #=======================================================================
        # # Sort results by ascending distance for each point
        # for coordinate in bounding_box: 
        #     point_result_dict[coordinate] = sorted(point_result_dict[coordinate], key=lambda d: d['distance'], reverse=False)       
        #=======================================================================
         
        logger.debug('Elapsed time: {}'.format(datetime.now() - t0))
    
def JW_metadata_filter(metadata_dict):
    '''
    Example function to filter datasets based on metadata values in metadata_dict
    This version applies John Wilford's filter conditions
    Returns True for match, False otherwise
    '''
    try:
        # Reject any datasets earlier than 1981
        if datetime.strptime(metadata_dict['acquisition_start_date'], '%d/%m/%Y') < datetime.strptime('01/01/1981', '%d/%m/%Y'):
            return False
            
        # Only accept GA/BMR/AGSO datasets between 1981 and 1992
        if (datetime.strptime(metadata_dict['acquisition_start_date'], '%d/%m/%Y') < datetime.strptime('01/01/1992', '%d/%m/%Y')
            and metadata_dict['client'] not in ['Geoscience Australia',
                                                'BMR',
                                                'AGSO',
                                                'GA'
                                                ]
            ):
                return False
    except ValueError:
        logger.warning('WARNING: Unhandled date format: {}'.format(metadata_dict['acquisition_start_date']))
        return False 
    
    return True        
                

def main():
    '''
    main routine for quick and dirty testing
    '''
    def get_args():
        """
        Handles all the arguments that are passed into the script

        :return: Returns a parsed version of the arguments.
        """
        parser = argparse.ArgumentParser(description='Bulk update records in eCat')
        parser.add_argument("-b", "--bounding_box",
                            help='bounding box in dataset native CRS expressed as "<min_xord>,<min_yord>,<max_xord>,<max_yord>"',
                            dest="bounding_box",
                            required=True)
        parser.add_argument("-m", "--metadata_path",
                            help="Path to CSV file containing netCDF metadata (default is {})".format(GeophysPointFetcher.DEFAULT_METADATA_CSV_PATH),
                            dest="metadata_path",
                            default=GeophysPointFetcher.DEFAULT_METADATA_CSV_PATH)
        parser.add_argument("-k", "--keywords",
                            help="keywords in the form <keyword1>[,<keyword2>,...]",
                            dest="keywords",
                            required=False)
        parser.add_argument("-v", "--variable_names",
                            help="variable_names in the form <variable_names1>[,<variable_names2>,...]",
                            dest="variable_names",
                            required=False)        
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='Output debug information. Default is no debug info')

        return parser.parse_args()    
    
    def parse_bounding_box_string(bounding_box_string):
        '''
        Helper function to parse bounding box coordinates from string of form "<min_xord>,<min_yord>,<max_xord>,<max_yord>"
        '''
        try:
            bounding_box = [float(value) for value in re.sub('\s', '', bounding_box_string).split(',')]
            
            assert len(bounding_box) == 4
            
            assert bounding_box[0] < bounding_box[2] and bounding_box[1] < bounding_box[3]
            
            return bounding_box
        except:
            raise Exception('Invalid bounding box string "{}"'.format(bounding_box_string))
    
    # Start of main function    
    # Parse the arguments passed in on the command line
    args = get_args()

    # Set Logging level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        
    bounding_box = parse_bounding_box_string(args.bounding_box)
    logger.debug('bounding_box = {}'.format(bounding_box))
        
    keywords = [keyword.strip() for keyword in args.keywords.split(',')] if args.keywords else []
    logger.debug('keywords = {}'.format(keywords))

    variable_names = [variable_name.strip() for variable_name in args.variable_names.split(',')] if args.variable_names else None
    logger.debug('variable_names = {}'.format(variable_names))

    gpf = GeophysPointFetcher()
    #logger.debug(gpf.__dict__)
    
    for dataset_dict in gpf.point_data_generator(bounding_box=bounding_box,
                                                 keywords=keywords,
                                                 metadata_filter_function=JW_metadata_filter,
                                                 variable_names=variable_names
                                                 ):
    
        logger.debug(pformat(dataset_dict))
    #logger.info('A total of {} points were found near {} in {} datasets'.format(len(point_list), coordinate, len(dataset_set)))

   

if __name__ == '__main__':
    main()
    
