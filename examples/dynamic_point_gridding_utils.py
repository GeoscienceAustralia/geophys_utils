'''
dynamic_point_gridding_utils.py
Implements utilities used for dynamic point dataset search and gridding 

Created on 3 May 2019

@author: Alex Ip
'''

import numpy as np
from netCDF4 import Dataset
import math
from scipy.interpolate import griddata
from geophys_utils import CSWUtils
from geophys_utils import NetCDFPointUtils
from geophys_utils import transform_coords, date_string2datetime

def get_netcdf_datasets(keywords, 
                        bounding_box=None, 
                        start_date_string=None, 
                        end_date_string=None, 
                        csw_url=None,
                        ):
    '''
    Find all datasets of interest and return a list of NetCDF file paths or OPeNDAP web service endpoints
    '''    
    csw_url = csw_url or 'https://ecat.ga.gov.au/geonetwork/srv/eng/csw'
    #create a csw_utils object and populate the parameters with search parameters
    try:
        cswu = CSWUtils(csw_url) 
    except:
        cswu = CSWUtils(csw_url, verify=False) 
        
    if start_date_string:
        start_datetime = date_string2datetime(start_date_string)
        assert start_datetime is not None, 'Invalid date string for start date'
    else:
        start_datetime = None
        
    if end_date_string:
        end_datetime = date_string2datetime(end_date_string)
        assert end_datetime is not None, 'Invalid date string for end date'
    else:
        end_datetime = None
        
    print('Querying CSW')
    record_list = [record for record in cswu.query_csw(keyword_list=keywords,
                                      #anytext_list=allwords,
                                      #titleword_list=titlewords,
                                      bounding_box=bounding_box,
                                      start_datetime=start_datetime,
                                      stop_datetime=end_datetime,
                                      #max_total_records=2000
                                      )
              ]
    print('{} matching dataset records found from CSW'.format(len(record_list)))
    
    netcdf_list = [distribution['url']
            for distribution in cswu.get_netcdf_urls(record_list)
            ]

    print('{} NetCDF distributions found'.format(len(netcdf_list)))
    
    return netcdf_list


def dataset_value_generator(variable_name_list,
                            dataset_list, 
                            bounding_box, 
                            min_points=None,
                            max_points=None
                            ):
    '''
    Generator yielding coordinates and values of the specified variable for all points from the supplied dataset list 
    which fall within bounds
    '''    
    line_dataset_count = 0
    for dataset in dataset_list:
        line_data = {}
        try:
            nc_dataset = Dataset(dataset + '#fillmismatch') # Note work-around for bad _FillValue: https://github.com/Unidata/netcdf-c/issues/1299
            netcdf_point_utils = NetCDFPointUtils(nc_dataset) 
            
            #print netcdf_point_utils.__dict__
            #print(nc_dataset.variables.keys())
            
            #print('Computing spatial mask')
            spatial_mask = netcdf_point_utils.get_spatial_mask(bounding_box)
            
            point_count = np.count_nonzero(spatial_mask)

            print('{}/{} points found in bounding box for {}'.format(point_count, netcdf_point_utils.point_count, dataset))
                        
            if not point_count:
                continue
            
            # Enforce min/max point counts
            if min_points and point_count < min_points:
                print('Skipping dataset with < {} points'.format(min_points))
                continue
            if max_points and point_count > max_points:
                print('Skipping dataset with > {} points'.format(max_points))
                continue
                
            dataset_value_dict = {'coordinates': netcdf_point_utils.xycoords[spatial_mask]}
            
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


def grid_points(coordinates,
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
    pixel_centre_bounds = (round(math.floor(grid_bounds[0] / grid_resolution) * grid_resolution, 6),
                   round(math.floor(grid_bounds[1] / grid_resolution) * grid_resolution, 6),
                   round(math.floor(grid_bounds[2] / grid_resolution - 1.0) * grid_resolution + grid_resolution, 6),
                   round(math.floor(grid_bounds[3] / grid_resolution - 1.0) * grid_resolution + grid_resolution, 6)
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

#===============================================================================
#     print("Computing spatial mask")
#     spatial_subset_mask = np.logical_and(np.logical_and((grid_bounds[0] <= grid_coordinates[:,0]), 
#                                                         (grid_coordinates[:,0] <= grid_bounds[2])), 
#                                          np.logical_and((grid_bounds[1] <= grid_coordinates[:,1]), 
#                                                         (grid_coordinates[:,1] <= grid_bounds[3]))
#                                         )    
#     # Skip points to reduce memory requirements
#     print("Generating point subset mask")
#     point_subset_mask = np.zeros(shape=values.shape, dtype=bool)
#     point_subset_mask[0:-1:point_step] = True
#     point_subset_mask = np.logical_and(spatial_subset_mask, point_subset_mask)
#     assert point_subset_mask.any(), 'No points found within grid bounds %s' % grid_bounds
#     
#     grid_coordinates = grid_coordinates[point_subset_mask]
#===============================================================================

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


def hillshade(array, azimuth, angle_altitude, vertical_exaggeration=1):
    '''
    Hillshade function adapted from:
    http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html
    '''       
    x, y = np.gradient(array*vertical_exaggeration)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth * np.pi / 180.
    altituderad = angle_altitude * np.pi / 180.
 
    shaded = np.sin(altituderad) * np.sin(slope)\
     + np.cos(altituderad) * np.cos(slope)\
     * np.cos(azimuthrad - aspect)
    
    return 255 * (shaded + 1) / 2    