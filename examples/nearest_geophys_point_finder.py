'''
Created on 11 Jan. 2018

@author: Alex Ip

Quick-and-dirty point finding demo
'''
import os
import sys
import csv
import netCDF4
import argparse
import logging
import re
from pprint import pformat

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

class NearestGeophysPointFinder(object):
    '''
    NearestGeophysPointFinder class definition
    '''
    DEFAULT_METADATA_CSV_PATH = os.path.join(os.path.dirname(__file__), 'geophysics_line_nc_metadata.csv')
    OPENDAP_PATH_MAP = ('/g/data2/uc0', 'http://dapds00.nci.org.au/thredds/dodsC/uc0')

    def __init__(self, metadata_csv_path=None):
        '''
        Constructor
        '''
        metadata_csv_path = metadata_csv_path or NearestGeophysPointFinder.DEFAULT_METADATA_CSV_PATH
        
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
        
        
    def dataset_metadata_generator(self, coordinate_list, max_distance=None, metadata_filter_dict={}):
        '''
        Generator returning all metadata dicts near given coordinate
        '''
        #TODO: Replace this bodgy CSV-based code with a CSW catalogue based solution
        max_distance = max_distance or 0
        
        for metadata_dict in self._metadata:
            for coordinate in coordinate_list:
                if not (float(metadata_dict['geospatial_lon_min']) <= (coordinate[0] + max_distance)
                    and float(metadata_dict['geospatial_lat_min']) <= (coordinate[1] + max_distance)
                    and float(metadata_dict['geospatial_lon_max']) >= (coordinate[0] - max_distance)
                    and float(metadata_dict['geospatial_lat_max']) >= (coordinate[1] - max_distance)
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
                    break # Only return one record for any number of coordinates


    def get_nearest_point_data(self, coordinate_list, points_required=1, max_distance=None, metadata_filter_dict={}):
        '''
        Function returning list of nearest points closest to coordinate with attributes
        '''
        point_result_dict = {coordinate: []
                             for coordinate in coordinate_list
                             }
        for metadata_dict in self.dataset_metadata_generator(coordinate_list, 
                                                            max_distance=max_distance,
                                                            metadata_filter_dict=metadata_filter_dict
                                                            ):
            nc_path = metadata_dict['file_path']
            
            if not os.path.isfile(nc_path): 
                nc_path = nc_path.replace(NearestGeophysPointFinder.OPENDAP_PATH_MAP[0], 
                                          NearestGeophysPointFinder.OPENDAP_PATH_MAP[1]
                                          )
                
            try:
                logger.info('Opening {}'.format(nc_path))
                nc_dataset = netCDF4.Dataset(nc_path, 'r')
                netcdf_line_utils = NetCDFLineUtils(nc_dataset)
                
                for coordinate in coordinate_list:
                    point_result_list = point_result_dict[coordinate]
                    distances, point_indices = netcdf_line_utils.nearest_neighbours(coordinate, points_required=points_required, max_distance=max_distance)
                    #logger.debug('distances = {}'.format(distances))
                    #logger.debug('point_indices = {}'.format(point_indices))
                    
                    if point_indices is not None:
                        # Convert scalars to lists if required
                        if not hasattr(point_indices, '__iter__'):
                            distances = [distances]
                            point_indices = [point_indices]
                            
                        logger.info('{} points near {} found in {}'.format(len(point_indices), coordinate, nc_path))
                        for found_index in range(len(point_indices)):
                            point_metadata_dict = dict(metadata_dict)
                            point_metadata_dict['distance'] = distances[found_index]
                            point_metadata_dict['coordinate'] = tuple(netcdf_line_utils.xycoords[point_indices[found_index]])
                            
                            for variable_name in netcdf_line_utils.point_variables:
                                point_metadata_dict[variable_name] = nc_dataset.variables[variable_name][point_indices[found_index]]
                        
                            #logger.debug('\tpoint_metadata_dict = {}'.format(point_metadata_dict))
                            point_result_list.append(point_metadata_dict)  
                    else:  
                        logger.info('No points near {} found in {}'.format(coordinate, nc_path))
            except Exception as e:
                logger.error('Failed to retrieve point data from {}: {}'.format(nc_path, e.message))
                
        # Sort results by ascending distance for each point
        for coordinate in coordinate_list: 
            point_result_dict[coordinate] = sorted(point_result_dict[coordinate], key=lambda d: d['distance'], reverse=False)       
        
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
        parser.add_argument("-c", "--coordinates",
                            help='coordinate pairs in dataset native CRS expressed as "(<xord>,<yord>)[,(<xord>,<yord>)...]"',
                            dest="coordinates",
                            required=True)
        parser.add_argument("-p", "--points_per_dataset",
                            help="Maximum number of points to retrieve from each dataset (default is 1)",
                            dest="points_per_dataset",
                            type=int,
                            default=1)
        parser.add_argument("-x", "--max_distance",
                            help="Maximum distance to search expressed in dataset native units (e.g. degrees, default is 0.1)",
                            dest="max_distance",
                            type=float,
                            default=0.1)
        parser.add_argument("-m", "--metadata_path",
                            help="Path to CSV file containing netCDF metadata (default is {})".format(NearestGeophysPointFinder.DEFAULT_METADATA_CSV_PATH),
                            dest="metadata_path",
                            default=NearestGeophysPointFinder.DEFAULT_METADATA_CSV_PATH)
        parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                            help='Output debug information. Default is no debug info')
        parser.add_argument('filter_args', nargs=argparse.REMAINDER,
                            help='Multiple filter arguments in the form of <key>=<value>')

        return parser.parse_args()    
    
    def parse_coordinate_list_string(coordinate_list_string):
        '''
        Helper function to parse coordinate list from string of form "(<xord>,<yord>)[,(<xord>,<yord>)...]"
        '''
        coordinate_list = []
        try:
            processing_string = coordinate_list_string.strip()
            while processing_string:
                match = re.match('\((([\+\-.0-9]+),(\s*)([\+\-.0-9]+))\)(,(.*))*', processing_string)
                if match:
                    coordinate_list.append((float(match.group(2)), float(match.group(4))))
                    if match.group(6) is not None:
                        processing_string = match.group(6).strip()
                    else:
                        processing_string = None
                    
            assert coordinate_list
            
            return coordinate_list
        except:
            raise Exception('Invalid coordinate list string "{}"'.format(coordinate_list_string))
    
    # Start of main function    
    # Parse the arguments passed in on the command line
    args = get_args()

    # Set Logging level
    if args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        
    coordinate_list = parse_coordinate_list_string(args.coordinates)

    logger.debug('coordinate_list = {}'.format(coordinate_list))
        
    logger.debug('points_per_dataset = {}'.format(args.points_per_dataset))
    logger.debug('max_distance = {}'.format(args.max_distance))

    # Set up modifier_args dict of key=value pairs
    metadata_filter_dict = {}
    for filter_arg in args.filter_args:
        match = re.match('^(\w+)=(.*)$', filter_arg)
        metadata_filter_dict[match.group(1)] = match.group(2)
    
    logger.debug('metadata_filter_dict = {}'.format(metadata_filter_dict))

    ngpf = NearestGeophysPointFinder()
    #logger.debug(ngpf.__dict__)
    
    point_result_dict = ngpf.get_nearest_point_data(coordinate_list=coordinate_list,
                                             points_required=args.points_per_dataset, 
                                             max_distance=args.max_distance,
                                             metadata_filter_dict=metadata_filter_dict
                                             )
    
    for coordinate, point_list in point_result_dict.iteritems():
        dataset_set = set()
        for point_dict in point_list:
            dataset_set.add(point_dict['file_path'])
            
        
        logger.info('A total of {} points were found near {} in {} datasets'.format(len(point_list), coordinate, len(dataset_set)))
        logger.debug(pformat(point_list[0:args.points_per_dataset]))

    
if __name__ == '__main__':
    main()
    
