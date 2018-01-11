'''
Created on 11 Jan. 2018

@author: u76345
'''
import os
import csv
import netCDF4
import numpy as np
from pprint import pprint

from geophys_utils import NetCDFLineUtils

class NearestGeophysPointFinder(object):
    '''
    classdocs
    '''
    DEFAULT_METADATA_CSV_PATH = os.path.join(os.path.dirname(__file__), 'geophysics_line_nc_metadata.csv')
    OPENDAP_PATH_MAP = ('/g/data2/uc0', 'http://dapds00.nci.org.au/thredds/dodsC/uc0')

    def __init__(self, metadata_csv_path=None):
        '''
        Constructor
        '''
        metadata_csv_path = metadata_csv_path or NearestGeophysPointFinder.DEFAULT_METADATA_CSV_PATH
        
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
        
        
    def metadata_containing_point(self, coordinates, max_distance=None):
        '''
        Generator returning all metadata dicts near given coordinates
        '''
        max_distance = max_distance or 0
        
        for metadata_dict in self._metadata:
            if (float(metadata_dict['geospatial_lon_min']) <= (coordinates[0] + max_distance)
                and float(metadata_dict['geospatial_lat_min']) <= (coordinates[1] + max_distance)
                and float(metadata_dict['geospatial_lon_max']) >= (coordinates[0] - max_distance)
                and float(metadata_dict['geospatial_lat_max']) >= (coordinates[1] - max_distance)
                ):
                yield metadata_dict


    def get_nearest_point_data(self, coordinates, points_required=1, max_distance=None):
        '''
        Function returning list of nearest points closest to coordinates with attributes
        '''
        point_info_list = []
        for metadata_dict in self.metadata_containing_point(coordinates, max_distance=max_distance):
            nc_path = metadata_dict['file_path']
            
            if not os.path.isfile(nc_path): 
                nc_path = nc_path.replace(NearestGeophysPointFinder.OPENDAP_PATH_MAP[0], 
                                          NearestGeophysPointFinder.OPENDAP_PATH_MAP[1]
                                          )
                
            print 'Opening {}'.format(nc_path)
            nc_dataset = netCDF4.Dataset(nc_path, 'r')
            netcdf_line_utils = NetCDFLineUtils(nc_dataset)
            
            distances, point_indices = netcdf_line_utils.nearest_neighbours(coordinates, points_required=points_required, max_distance=max_distance)
            if point_indices is not None:
                # Convert scalars to lists if required
                if not hasattr(point_indices, '__iter__'):
                    distances = [distances]
                    point_indices = [point_indices]
                    
                print '{} points found in {}'.format(len(point_indices), nc_path)
                for index in range(len(point_indices)):
                    point_metadata_dict = dict(metadata_dict)
                    point_metadata_dict['distance'] = distances[index]
                    point_metadata_dict['coordinates'] = netcdf_line_utils.xycoords[point_indices[index]]
                    
                    for variable_name in netcdf_line_utils.point_variables:
                        point_metadata_dict[variable_name] = nc_dataset.variables[variable_name][point_indices[index]]
                
                point_info_list.append(point_metadata_dict)  
            else:  
                print 'No points found in {}'.format(nc_path)
                
        return sorted(point_info_list, key=lambda d: d['distance'], reverse=False) 
    

def main():
    '''
    main routine for quick and dirty testing
    '''
    ngpf = NearestGeophysPointFinder()
    #pprint(ngpf.__dict__)
    pprint(ngpf.get_nearest_point_data(coordinates=(149, -35),
                                       points_required=10, 
                                       max_distance=0.2 
                                       )
           )

    
if __name__ == '__main__':
    main()
    
