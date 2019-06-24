'''
Created on 14 Jun 2019

@author: Alex Ip
'''
import numpy as np
import math
from geophys_utils import get_spatial_ref_from_wkt, get_wkt_from_spatial_ref, transform_coords
from netCDF4 import Dataset
from scipy.interpolate import griddata
from geophys_utils import CSWUtils
from geophys_utils import NetCDFPointUtils
from geophys_utils import array2file
import os
import sys
import re
import argparse
import itertools
from pprint import pprint

DEBUG = True

class TileGridder(object):
    '''
    TileGridder
    '''
    DEFAULT_TILE_EXPANSION_PERCENT = 0.05
    DEFAULT_CSW_URL = 'https://ecat.ga.gov.au/geonetwork/srv/eng/csw'
    GDA94_CRS_WKT = get_wkt_from_spatial_ref(get_spatial_ref_from_wkt('EPSG:4283')) # Defaults to GDA94
    DEFAULT_FILTER_VARIABLE_NAME = 'gridflag'
    DEFAULT_FILTER_VALUE_LIST = ['Station used in the production of GA grids.']

    def __init__(self, 
                 dataset_keywords, 
                 grid_variable_name, # Name of data variable to grid
                 grid_bounds, # Grid bounds as [xmin, ymin, xmax, ymax]
                 grid_crs_wkt = None, # Defaults to GDA94 
                 start_datetime=None,
                 end_datetime=None,
                 filter_variable_name=None, # e.g. 'gridflag'
                 filter_value_list=None, # e.g. ['Station used in the production of GA grids.']
                 tile_extra=None, # Absolute extra per side. Defaults to 5% extra on each side
                 ):
        '''
        TileGridder Constructor
        '''
        self._dataset_list = None
        self._dataset_values = None
        
        self.dataset_keywords = dataset_keywords
        self.grid_variable_name = grid_variable_name
        self.grid_bounds = grid_bounds
        self.grid_crs_wkt = get_wkt_from_spatial_ref(get_spatial_ref_from_wkt(grid_crs_wkt or TileGridder.GDA94_CRS_WKT))
        self.filter_variable_name = filter_variable_name or TileGridder.DEFAULT_FILTER_VARIABLE_NAME
        self.filter_value_list = filter_value_list or TileGridder.DEFAULT_FILTER_VALUE_LIST
        self.start_datetime = start_datetime
        self.end_datetime = end_datetime
        
        
        if tile_extra is None: # Expand bounds by TileGridder.DEFAULT_TILE_EXPANSION_PERCENT percent each side
            self.expanded_grid_bounds = [grid_bounds[0] - (grid_bounds[2] - grid_bounds[0]) * TileGridder.DEFAULT_TILE_EXPANSION_PERCENT,
                                         grid_bounds[1] - (grid_bounds[3] - grid_bounds[1]) * TileGridder.DEFAULT_TILE_EXPANSION_PERCENT,
                                         grid_bounds[2] + (grid_bounds[2] - grid_bounds[0]) * TileGridder.DEFAULT_TILE_EXPANSION_PERCENT,
                                         grid_bounds[3] + (grid_bounds[3] - grid_bounds[1]) * TileGridder.DEFAULT_TILE_EXPANSION_PERCENT,
                                         ]
        else: # Expand bounds by absolute amount
            self.expanded_grid_bounds = [grid_bounds[0] - tile_extra,
                                         grid_bounds[1] - tile_extra,
                                         grid_bounds[2] + tile_extra,
                                         grid_bounds[3] + tile_extra,
                                         ]

        self.expanded_gda94_grid_bounds = self.reproject_bounds(self.expanded_grid_bounds,
                                                                self.grid_crs_wkt,
                                                                TileGridder.GDA94_CRS_WKT)
        
        print('expanded_gda94_grid_bounds = {}'.format(self.expanded_gda94_grid_bounds))


    @property
    def dataset_list(self):
        '''
        list of individual datasets which intersect bounding box
        '''
        if self._dataset_list is None:
        
            self._dataset_list = self.get_netcdf_datasets(self.dataset_keywords, 
                                                         bounding_box=self.expanded_gda94_grid_bounds, 
                                                         start_datetime=self.start_datetime, 
                                                         end_datetime=self.end_datetime, 
                                                         csw_url=None,
                                                         )
        
        return self._dataset_list

    @property
    def dataset_values(self):
        '''
        Read and filter points from individual datasets
        '''
        if self._dataset_values is None:
            self._dataset_values = {dataset: dataset_value_dict
                for dataset, dataset_value_dict in self.dataset_value_generator([self.grid_variable_name, self.filter_variable_name],
                                                                                self.dataset_list, 
                                                                                self.expanded_gda94_grid_bounds,
                                                                                )
                                    }
            self.filter_points()
            
        return self._dataset_values
        

    def filter_points(self):        
        '''
        Set filter points from individual datasets
        e.g. Only use points where gridflag == 'Station used in the production of GA grids.'
        '''
        # Only filter if we have a filter variable and allowed values
        if not (self.filter_variable_name and self.filter_value_list):
            return
        
        for dataset in sorted(self.dataset_values.keys()):
            filter_mask = np.zeros(shape=(self.dataset_values[dataset]['coordinates'].shape[0],), dtype=np.bool)
            for filter_value in self.filter_value_list:
                filter_mask = np.logical_or(filter_mask, (self.dataset_values[dataset][self.filter_variable_name] == filter_value))
                
            coordinates = self.dataset_values[dataset]['coordinates'][filter_mask]
            if len(coordinates):
                self.dataset_values[dataset]['coordinates'] = coordinates
                self.dataset_values[dataset][self.grid_variable_name] = self.dataset_values[dataset][self.grid_variable_name][filter_mask]
                del self.dataset_values[dataset][self.filter_variable_name] # We don't need this any more
            else:
                del self.dataset_values[dataset] # No usable points in dataset
                
            
    def reproject_bounds(self, bounds, from_crs_wkt, to_crs_wkt):
        '''
        Function to return orthogonal bounds reprojected to new CRS
        '''
        if from_crs_wkt == to_crs_wkt: # No change
            return bounds
        
        bounding_box = ((bounds[0], bounds[1]), 
                        (bounds[2], bounds[1]), 
                        (bounds[2], bounds[3]), 
                        (bounds[0], bounds[3])
                        )
        
        reprojected_bounding_box = np.array(transform_coords(bounding_box, 
                                                             from_crs_wkt, 
                                                             to_crs_wkt))
        
        reprojected_bounds = (min(reprojected_bounding_box[:,0]), 
                              min(reprojected_bounding_box[:,1]), 
                              max(reprojected_bounding_box[:,0]), 
                              max(reprojected_bounding_box[:,1])
                              )
        
        return reprojected_bounds
    

    def get_netcdf_datasets(self,
                            keywords, 
                            bounding_box=None, 
                            start_datetime=None, 
                            end_datetime=None, 
                            csw_url=None,
                            ):
        '''
        Find all datasets of interest and return a list of NetCDF file paths or OPeNDAP web service endpoints
        '''    
        csw_url = csw_url or TileGridder.DEFAULT_CSW_URL
        #create a csw_utils object and populate the parameters with search parameters
        try:
            cswu = CSWUtils(csw_url) 
        except:
            cswu = CSWUtils(csw_url, verify=False) 
            
        print('Querying CSW')
        record_list = [record for record in cswu.query_csw(keyword_list=keywords,
                                          #anytext_list=allwords,
                                          #titleword_list=titlewords,
                                          bounding_box=bounding_box,
                                          start_datetime=start_datetime,
                                          stop_datetime=end_datetime,
                                          #max_total_records=2000,
                                          get_layers=False,
                                          )
                  ]
        print('{} matching dataset records found from CSW'.format(len(record_list)))
        
        netcdf_list = [str(distribution['url'])
                for distribution in cswu.get_netcdf_urls(record_list)
                ]
    
        print('{} NetCDF distributions found'.format(len(netcdf_list)))
        
        return netcdf_list
    
    
    def dataset_value_generator(self,
                                variable_name_list,
                                dataset_list, 
                                bounding_box, 
                                min_points=None,
                                max_points=None
                                ):
        '''
        Generator yielding coordinates and values of the specified variable for all points from the supplied dataset list 
        which fall within bounds
        '''    
        for dataset in dataset_list:
            try:
                try:
                    nc_dataset = Dataset(dataset)
                except:
                    nc_dataset = Dataset(dataset + '#fillmismatch') # Note work-around for bad _FillValue: https://github.com/Unidata/netcdf-c/issues/1299
                
                netcdf_point_utils = NetCDFPointUtils(nc_dataset) 
                
                #print netcdf_point_utils.__dict__
                #print(nc_dataset.variables.keys())
                #print('Computing spatial mask')
                spatial_mask = netcdf_point_utils.get_spatial_mask(bounding_box, self.grid_crs_wkt)
                point_count = np.count_nonzero(spatial_mask)
    
                print('{}/{} points found in expanded bounding box for {}'.format(point_count, netcdf_point_utils.point_count, dataset))
                            
                if not point_count:
                    continue
                
                # Enforce min/max point counts
                if min_points and point_count < min_points:
                    print('Skipping dataset with < {} points'.format(min_points))
                    continue
                if max_points and point_count > max_points:
                    print('Skipping dataset with > {} points'.format(max_points))
                    continue
                    
                dataset_value_dict = {'coordinates': transform_coords(netcdf_point_utils.xycoords[spatial_mask],
                                                                      get_wkt_from_spatial_ref(get_spatial_ref_from_wkt(netcdf_point_utils.wkt)),
                                                                      self.grid_crs_wkt
                                                                      )
                                      }
                
                # Read all variable attributes and values
                for variable_name in variable_name_list:
                    variable = nc_dataset.variables[variable_name]
                    if (variable.dimensions[0] != 'point'): # Variable is NOT of point dimension - must be lookup
                        dataset_value_dict[variable_name] = netcdf_point_utils.expand_lookup_variable(lookup_variable_name=variable_name, 
                                                                                                      mask=spatial_mask)                     
                    else: # 'point' is in variable.dimensions - "normal" variable                
                        dataset_value_dict[variable_name] = variable[spatial_mask]
                
                yield dataset, dataset_value_dict
        
            except Exception as e:
                print('Unable to read point dataset {}: {}'.format(dataset, e))
    
    
    def grid_points(self,
                    coordinates,
                    coordinate_wkt,
                    values,
                    grid_wkt, 
                    grid_bounds,
                    grid_resolution, 
                    resampling_method='linear', 
                    point_step=1):
        '''
        Return geotransform CRS WKT, and interpolated grid from supplied coordinates and points
        '''
        
        # Determine spatial grid bounds rounded out to nearest GRID_RESOLUTION multiple
        pixel_centre_bounds = (round((math.floor(grid_bounds[0] / grid_resolution) + 0.5) * grid_resolution, 6),
                       round((math.floor(grid_bounds[1] / grid_resolution) + 0.5) * grid_resolution, 6),
                       round((math.floor(grid_bounds[2] / grid_resolution) - 0.5) * grid_resolution, 6),
                       round((math.floor(grid_bounds[3] / grid_resolution) - 0.5) * grid_resolution, 6)
                       )
        
        print("Reprojecting coordinates")
        grid_coordinates = np.array(transform_coords(coordinates, coordinate_wkt, grid_wkt))
        #print('grid_coordinates = {}'.format(grid_coordinates))
    
        # Create grids of Y and X values. Note YX ordering and inverted Y for image
        # Note GRID_RESOLUTION/2.0 fudge to avoid truncation due to rounding error
        print("Generating grid coordinates")
        grid_y, grid_x = np.mgrid[pixel_centre_bounds[3]:pixel_centre_bounds[1]-grid_resolution/2.0:-grid_resolution, 
                                  pixel_centre_bounds[0]:pixel_centre_bounds[2]+grid_resolution/2.0:grid_resolution]
     
        
        # Skip points to reduce memory requirements
        print("Generating point subset mask")
        point_subset_mask = np.zeros(shape=values.shape, dtype=bool)
        point_subset_mask[0:-1:point_step] = True
         
        grid_coordinates = grid_coordinates[point_subset_mask]
        values = values[point_subset_mask]
    
        # Interpolate required values to the grid - Note yx ordering and inverted y for image
        print("Interpolating {} points".format(grid_coordinates.shape[0]))
        grid_array = griddata(grid_coordinates[:,::-1],
                              values,
                              (grid_y, grid_x), 
                              method=resampling_method
                              )
    
        print("Interpolation complete")
        #  crs:GeoTransform = "109.1002342895272 0.00833333 0 -9.354948067227777 0 -0.00833333 "
        geotransform = [pixel_centre_bounds[0]-grid_resolution/2.0,
                        grid_resolution,
                        0,
                        pixel_centre_bounds[3]+grid_resolution/2.0,
                        0,
                        -grid_resolution
                        ] 
    
        return grid_array, grid_wkt, geotransform
        
        
    def grid_tile(self, 
                  grid_resolution, 
                  coordinates=None,
                  values=None,
                  resampling_method='linear', 
                  point_step=1):
        
        if coordinates is None:
            coordinates = np.concatenate([self.dataset_values[dataset]['coordinates']
                                          for dataset in sorted(self.dataset_values.keys())
                                          ]
                                         )
                                         
        if values is None:
            values = np.concatenate([self.dataset_values[dataset][self.grid_variable_name]
                                          for dataset in sorted(self.dataset_values.keys())
                                          ]
                                         )
                                         
        return self.grid_points(coordinates=coordinates,
                                coordinate_wkt=self.grid_crs_wkt,
                                values=values,
                                grid_wkt=self.grid_crs_wkt, 
                                grid_bounds=self.grid_bounds,
                                grid_resolution=grid_resolution, 
                                resampling_method=resampling_method, 
                                point_step=point_step)
        
        
    def output_points(self, point_list_path):
        '''
        Write CSV containing all points to point_list_path
        '''
        with open(point_list_path, 'w') as output_file:
            output_file.write(', '.join(['Dataset', 'X', 'Y', self.grid_variable_name]) + '\n')
            for dataset in sorted(self.dataset_values.keys()):
                #for array_name in sorted(self.dataset_values[dataset].keys()):
                arrays = self.dataset_values[dataset]
                for point_index in range(len(arrays['coordinates'])):
                    output_file.write(','.join([os.path.basename(dataset),
                                                    ','.join([str(ordinate) for ordinate in arrays['coordinates'][point_index]]),
                                                    str(arrays[self.grid_variable_name][point_index])
                                                    ]
                                                    ) + '\n'
                                       )
                                             
                
    def output_dataset_list(self, dataset_list_path):
        '''
        Write a text file containing all dataset paths or URLs to dataset_list_path
        '''
        with open(dataset_list_path, 'w') as output_file:
            for dataset in sorted(self.dataset_list):
                output_file.write(dataset + '\n')
            
    def read_dataset_list(self, dataset_list_path):
        '''
        Write a text file containing all dataset paths or URLs to dataset_list_path
        '''
        with open(dataset_list_path, 'r') as input_file:
            self._dataset_list = list([dataset.strip() for dataset in input_file.readlines()])
            
        
def main():
    '''
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
            
        
    #===========================================================================
    # dataset_keywords = 'point, gravity, point located data, ground digital data, geophysical survey' 
    # grid_variable_name = 'bouguer' # Name of data variable to grid
    # #grid_bounds = (118.75, -28.5, 119.5, -27.75) # Sandstone, WA
    # grid_bounds = (141.0, -32.5, 142.0, -31.5) # Broken Hill, NSW (~40k points - takes a while)
    # #grid_bounds = (136.5, -31.0, 137.5, -30.0) # Roxby Downs, SA (possibly some levelling issues?)
    # grid_crs_wkt = get_wkt_from_spatial_ref(get_spatial_ref_from_wkt('EPSG:4283')) # Defaults to GDA94 
    # start_datetime = None
    # end_datetime = None
    # filter_variable_name = 'gridflag'
    # filter_value_list = ['Station used in the production of GA grids.']
    # tile_extra=None # Absolute extra per side. Defaults to 5% extra on each side
    # grid_resolution=0.001
    # resampling_method='cubic'
    #===========================================================================

    # Define command line arguments
    parser = argparse.ArgumentParser()
    # Required arguments
    parser.add_argument("-k", "--keywords", help="comma-separated list of required keywords for search", type=str, required=True)
    parser.add_argument("-b", "--bounds", help='comma-separated <minx>,<miny>,<maxx>,<maxy> ordinates of bounding box for gridding. N.B: A leading "-" sign on this list should NOT be preceded by a space',
                        type=str, required=True)
    parser.add_argument("-t", "--tilesize", help="single value or comma-separated pair of x,y values for tile size", type=str, required=True)
    parser.add_argument("-d", "--data_variable", help="data variable name", type=str, required=True)
    parser.add_argument("-r", "--grid_resolution", help='grid resolution', type=float, required=True)
    
    parser.add_argument("-c", "--crs", help='coordinate reference system for bounding box coordinates for search. Defaults to "EPSG:4283".',
                        type=str, default='EPSG:4283')
    parser.add_argument("-f", "--filter_variable", help='name of filter variable. Defaults to "gridflag"', type=str, default='gridflag')    
    parser.add_argument("-v", "--filter_values", 
                        help='comma separated list of allowed filter values. Defaults to "Station used in the production of GA grids."', 
                        type=str, default='Station used in the production of GA grids.')
    parser.add_argument("-m", "--resampling_method", help='resampling method. Defaults to "cubic"', type=str)
    parser.add_argument("-x", "--tile_extra", help='absolute extra per side for tiling. Defaults to 5% extra on each side', type=float, 
                        required=False)
    #parser.add_argument("-s", "--start_date", help="start date for search", type=str)
    #parser.add_argument("-e", "--end_date", help="end date for search", type=str)
    parser.add_argument('-p', '--process', action='store_const', const=True, default=False,
                        help='Process tiles. Default is not to process, but just to create point CSV files')
    parser.add_argument('--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("-o", "--output_dir", help='output directory. Defaults to ".".', type=str, default='.')   
    args = parser.parse_args()


    dataset_keywords = args.keywords
    grid_variable_name = args.data_variable
    grid_resolution = args.grid_resolution
    grid_bounds = [round(float(ordinate.strip()) / grid_resolution) * grid_resolution for ordinate in args.bounds.split(',')]
    tile_size = [round(float(size.strip()) / grid_resolution) * grid_resolution for size in args.tilesize.split(',')]
    if len(tile_size) == 1: # Only one size given
        tile_size.append(tile_size[0]) # Use same size for X and Y
    
    grid_crs_wkt = get_wkt_from_spatial_ref(get_spatial_ref_from_wkt(args.crs)) # Defaults to GDA94
    start_datetime = None # args.start_date
    end_datetime = None # args.end_date
    filter_variable_name = args.filter_variable
    filter_value_list = [value.strip() for value in args.filter_values.split(',')]
    tile_extra = args.tile_extra # Absolute extra per side. Defaults to 5% extra on each side
    resampling_method = os.environ.get('filter_variable_name') or 'cubic'
        
    assert os.path.isdir(args.output_dir), 'Invalid output directory'
    output_dir = os.path.abspath(args.output_dir)
    
    os.makedirs(os.path.join(output_dir, 'tiles'))
    os.makedirs(os.path.join(output_dir, 'point_lists'))
    os.makedirs(os.path.join(output_dir, 'dataset_lists'))
    
    for ll_point in itertools.product(*[np.arange(grid_bounds[0+dim_index], grid_bounds[2+dim_index], tile_size[dim_index]) for dim_index in range(2)]):
        tile_bounds = [ll_point[0], ll_point[1], ll_point[0]+tile_size[0], ll_point[1]+tile_size[1]]
        print(tile_bounds)
    
        tile_path = os.path.join(output_dir, 'tiles', '_'.join([grid_variable_name] + [str(ordinate) for ordinate in tile_bounds]) + '.tif')
        point_list_path = os.path.join(output_dir, 'point_lists', '_'.join([grid_variable_name] + [str(ordinate) for ordinate in tile_bounds]) + '.csv')
        dataset_list_path = os.path.join(output_dir, 'dataset_lists', '_'.join([grid_variable_name] + [str(ordinate) for ordinate in tile_bounds]) + '.txt')
        
        tg = TileGridder(dataset_keywords, 
                         grid_variable_name, # Name of data variable to grid
                         tile_bounds, # Grid bounds as [xmin, ymin, xmax, ymax]
                         grid_crs_wkt, # Defaults to GDA94 
                         start_datetime,
                         end_datetime,
                         filter_variable_name, # e.g. 'gridflag'
                         filter_value_list, # e.g. ['Station used in the production of GA grids.']
                         tile_extra, # Absolute extra per side. Defaults to 5% extra on each side
                         )
        
        
        if not os.path.isfile(dataset_list_path): # No dataset list file exists - try to make one
            try:
                # pprint(tg.dataset_values)
                tg.output_dataset_list(dataset_list_path)
                print('Finished writing dataset list file {}'.format(dataset_list_path))
            except Exception as e:
                print('Unable to create dataset list file {}: {}'.format(dataset_list_path, e))
                
        # Process tiles if required
        if args.process and os.path.isfile(dataset_list_path):
            # Skip processing tile if already done
            if os.path.isfile(tile_path):
                print('Already created tile file {}'.format(tile_path))
                continue
            
            if os.path.isfile(dataset_list_path): # No point file exists - try to make one
                try:
                    # pprint(tg.dataset_values)
                    tg.read_dataset_list(dataset_list_path)
                    print('Finished read dataset list file {}'.format(dataset_list_path))
                except Exception as e:
                    print('Unable to read dataset list file {}: {}'.format(dataset_list_path, e))
                    continue
                
            if not len(tg.dataset_list):
                print('No datasets to process')
                continue
            
            if not os.path.isfile(point_list_path): # No point file exists - try to make one
                try:
                    # pprint(tg.dataset_values)
                    tg.output_points(point_list_path)
                    print('Finished writing point file {}'.format(point_list_path))
                except Exception as e:
                    print('Unable to create point file {}: {}'.format(point_list_path, e))
                    continue
                
            # Process tile
            point_coordinates = []
            point_values = []
            with open(point_list_path, 'r') as csv_file:
                _header = csv_file.readline() # Read header
                
                for line in csv_file.readlines():
                    line_values = [value.strip() for value in line.split(',')]
                    point_coordinates.append(tuple(float(ordinate) for ordinate in line_values[1:3]))
                    point_values.append(float(line_values[3]))
                                  
            point_coordinates = np.array(point_coordinates)
            point_values = np.array(point_values)
            print('Finished reading point file {}'.format(point_list_path))
            print(point_coordinates)
                    
            if not len(point_values):
                print('No points to grid')
                continue
            
            grid_array, grid_wkt, geotransform = tg.grid_tile(grid_resolution=grid_resolution, 
                                                              coordinates=point_coordinates,
                                                              values=point_values,
                                                              resampling_method=resampling_method, 
                                                              point_step=1)
             
            print(grid_array.shape, grid_wkt, geotransform)
             
            array2file(data_arrays=[grid_array], 
                       projection=grid_crs_wkt, 
                       geotransform=geotransform, 
                       file_path=tile_path, 
                       file_format='GTiff')

            

    
    
if __name__ == '__main__':
    main()
