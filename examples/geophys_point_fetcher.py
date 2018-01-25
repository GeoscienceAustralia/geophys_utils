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
        
        
    def dataset_metadata_generator(self, bounding_box, metadata_filter_dict={}, date_from=None):
        '''
        Generator returning all metadata dicts near given coordinate
        '''
        #TODO: Replace this bodgy CSV-based code with a CSW catalogue based solution
        for metadata_dict in self._metadata:
            if date_from:
                try:
                    if date_from > datetime.strptime(metadata_dict['acquisition_start_date'], '%d/%m/%Y'):
                        continue
                except ValueError:
                    logger.warning('WARNING: Unhandled date format: {}'.format(metadata_dict['acquisition_start_date']))
                    continue 
            
            if not (float(metadata_dict['geospatial_lon_min']) <= bounding_box[0]
                    and float(metadata_dict['geospatial_lat_min']) <= bounding_box[1]
                    and float(metadata_dict['geospatial_lon_max']) >= bounding_box[2]
                    and float(metadata_dict['geospatial_lat_max']) >= bounding_box[3]
                    ):
                continue
                
            # Yield only dataset metadata matching all key=value pairs in metadata_filter_dict
            match_found = True
            for key, value in metadata_filter_dict.iteritems():
                match_found = match_found and metadata_dict.get(key) == value
                if not match_found:
                    break
                
            if match_found:
                yield metadata_dict


    def get_point_data(self, bounding_box, metadata_filter_dict={}, date_from=None):
        '''
        Function returning list of nearest points closest to coordinate with attributes
        '''
        t0 = datetime.now()
        point_result_dict = {}
        for metadata_dict in self.dataset_metadata_generator(bounding_box, 
                                                            metadata_filter_dict=metadata_filter_dict,
                                                            date_from=date_from
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
                
                logger.debug('Computing spatial subset mask...')
                spatial_mask = netcdf_line_utils.get_spatial_mask(bounding_box)
                if not np.any(spatial_mask):
                    logger.warning('No points in bounding box {} for {}'.format(tuple(bounding_box), nc_path))
                    continue
                
                point_indices = np.where(spatial_mask)[0]
                logger.info('{} points found in bounding box {} for {}'.format(len(point_indices), tuple(bounding_box), nc_path))

                dataset_dict = {}
                
                dataset_dict['metadata'] = dict(metadata_dict)
                logger.debug('\tReading coordinates')
                dataset_dict['coordinates'] = netcdf_line_utils.xycoords[point_indices]
                for variable_name in netcdf_line_utils.point_variables:
                    logger.debug('\tReading values for {}'.format(variable_name))
                    dataset_dict[variable_name] = nc_dataset.variables[variable_name][point_indices]
                
                point_result_dict[nc_path] = dataset_dict
            except Exception as e:
                logger.error('Failed to retrieve point data from {}: {}'.format(nc_path, e.message))
                
        #=======================================================================
        # # Sort results by ascending distance for each point
        # for coordinate in bounding_box: 
        #     point_result_dict[coordinate] = sorted(point_result_dict[coordinate], key=lambda d: d['distance'], reverse=False)       
        #=======================================================================
         
        logger.debug('Elapsed time: {}'.format(datetime.now() - t0))
        return point_result_dict 
    

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
        parser.add_argument("-f", "--date_from",
                            help="Oldest date for survey in the form dd/mm/yyyy",
                            dest="date_from",
                            required=False)
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='Output debug information. Default is no debug info')
        parser.add_argument('filter_args', nargs=argparse.REMAINDER,
                            help='Multiple filter arguments in the form of <key>=<value>')

        return parser.parse_args()    
    
    def parse_bounding_box_string(bounding_box_string):
        '''
        Helper function to parse coordinate list from string of form "(<xord>,<yord>)[,(<xord>,<yord>)...]"
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
        
    date_from = None
    if args.date_from:
        date_from = datetime.strptime(args.date_from, '%d/%m/%Y')
        logger.debug('date_from = {}'.format(date_from.isoformat()))
        

    # Set up modifier_args dict of key=value pairs
    metadata_filter_dict = {}
    for filter_arg in args.filter_args:
        match = re.match('^(\w+)=(.*)$', filter_arg)
        metadata_filter_dict[match.group(1)] = match.group(2)
    
    logger.debug('metadata_filter_dict = {}'.format(metadata_filter_dict))

    gpf = GeophysPointFetcher()
    #logger.debug(gpf.__dict__)
    
    point_result_dict = gpf.get_point_data(bounding_box=bounding_box,
                                           metadata_filter_dict=metadata_filter_dict,
                                           date_from=date_from
                                           )
    
    #logger.info('A total of {} points were found near {} in {} datasets'.format(len(point_list), coordinate, len(dataset_set)))
    logger.debug(pformat(point_result_dict))

    
if __name__ == '__main__':
    main()
    
