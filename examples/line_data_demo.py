import math
import os

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata

from geophys_utils import CSWUtils
from geophys_utils import NetCDFLineUtils
from geophys_utils import date_string2datetime
from geophys_utils import transform_coords, get_spatial_ref_from_wkt, get_reprojected_bounds

# Setup proxy as required
GA_STAFF_WIFI = False

if GA_STAFF_WIFI:
    os.environ['http_proxy'] = 'http://proxy.inno.lan:3128'
    os.environ['https_proxy'] = 'http://proxy.inno.lan:3128'

# N.B: GA internal CSW addresses will need port forwarding to work from the NCI
# Also, dev.public.ecat.ga.gov.au requires a hack to csw_utils and owslib to overcome certificate problem
# DEFAULT_CSW_URL = 'https://dev.public.ecat.ga.gov.au/geonetwork/srv/eng/csw'  # GA's internally-facing development eCat
DEFAULT_CSW_URL = 'https://ecat.ga.gov.au/geonetwork/srv/eng/csw'  # GA's externally-facing eCat
csw_url = 'https://ecat.ga.gov.au/geonetwork/srv/eng/csw'  # GA's externally-facing eCat
# csw_url = 'https://internal.ecat.ga.gov.au/geonetwork/srv/eng/csw' # GA's internally-facing eCat
# csw_url = 'http://geonetworkrr2.nci.org.au/geonetwork/srv/eng/csw' # NCI GeoNetwork

WGS84_WKT = get_spatial_ref_from_wkt('EPSG:4326').ExportToWkt()


# Set up search criteria
# bounds = (120.0, -29.0, 121, -28) # Spatial subset of dataset in WGS84 coordinates
# keywords = 'geophysics,airborne digital data,geophysical survey,magnetics,line,AWAGS' # Comma-separated list of keywords


# Set spatial information about bounds
# centre_coords = [(bounds[dim_index] + bounds[dim_index+2]) / 2.0 for dim_index in range(2)]

# utm_wkt = get_utm_wkt(centre_coords, wgs84_wkt)
# reprojected_bounding_box = np.array(transform_coords(((bounds[0], bounds[1]), (bounds[2], bounds[1]), (bounds[2], bounds[3]), (bounds[0], bounds[3])), wgs84_wkt, utm_wkt))
# utm_bounds = [min(reprojected_bounding_box[:,0]),
#              min(reprojected_bounding_box[:,1]),
#              max(reprojected_bounding_box[:,0]),
#              max(reprojected_bounding_box[:,1])]

# print wgs84_wkt
# print centre_coords
# print utm_wkt
# print utm_bounds


def get_netcdf_datasets(keywords,
                        wgs84_bounds=None,
                        start_date_string=None,
                        end_date_string=None,
                        csw_url=None):
    '''
    Find all datasets of interest and return a list of NetCDF file paths or OPeNDAP web service endpoints
    '''
    csw_url = csw_url or DEFAULT_CSW_URL
    # create a csw_utils object and populate the parameters with search parameters
    cswu = CSWUtils(csw_url)

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
                                                       # anytext_list=allwords,
                                                       # titleword_list=titlewords,
                                                       bounding_box=wgs84_bounds,
                                                       start_datetime=start_datetime,
                                                       stop_datetime=end_datetime,
                                                       # max_total_records=2000
                                                       )
                   ]
    print(f'{len(record_list)} matching dataset records found from CSW')

    netcdf_list = [distribution['url']
                   for distribution in cswu.get_netcdf_urls(record_list)
                   ]

    print(f'{len(netcdf_list)} NetCDF distributions found')

    return netcdf_list


def dataset_point_generator(dataset_list,
                            wgs84_bounds,
                            variable_name,
                            coordinate_wkt,
                            flight_lines_only=True,
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
        print(f'Reading and reprojecting points from line dataset {dataset}')
        try:
            nc_dataset = Dataset(dataset)
            mag_awags_variable = nc_dataset.variables[variable_name]
            netcdf_line_utils = NetCDFLineUtils(nc_dataset)

            reprojected_bounds = get_reprojected_bounds(wgs84_bounds, WGS84_WKT,
                                                        netcdf_line_utils.wkt)
            # print netcdf_line_utils.__dict__

            if flight_lines_only:
                print('Excluding tie-lines')
                line_numbers = nc_dataset.variables['line'][nc_dataset.variables['lineType'][:] == 2]
                line_mask = np.zeros(shape=nc_dataset.variables[variable_name].shape, dtype=bool)
                for _line_number, single_line_mask in netcdf_line_utils.get_line_masks(line_numbers):
                    line_mask = np.logical_or(line_mask, single_line_mask)
            else:
                line_mask = np.ones(shape=nc_dataset.variables[variable_name].shape, dtype=bool)

            print('Computing spatial mask')
            selection_mask = np.logical_and(netcdf_line_utils.get_spatial_mask(reprojected_bounds),
                                            line_mask
                                            )
            selection_count = np.count_nonzero(selection_mask)
            print(f'{selection_count}/{len(mag_awags_variable)} points found in bounding box')

            # Enforce min/max point counts
            if min_points and selection_count < min_points:
                print(f'Skipping dataset with < {min_points} points')
                continue
            if max_points and selection_count > max_points:
                print(f'Skipping dataset with > {max_points} points')
                continue

            print('Transforming coordinates')
            coordinates = np.array(transform_coords(netcdf_line_utils.xycoords[selection_mask],
                                                    netcdf_line_utils.wkt,
                                                    coordinate_wkt))
            print(f'Transformed coordinates: {coordinates}')

            print('Extracting values')
            values = mag_awags_variable[selection_mask]
            print(f'Extracted values: {values}')

            print('Checking mask')
            mask = np.ma.getmask(values)
            if mask is not np.ma.nomask:
                print(f'Discarding {np.count_nonzero(mask)} invalid values')
                values = values[~mask].data
                coordinates = coordinates[~mask]
                print(f"{values.shape[0]} valid points were found")

            line_dataset_count += 1
            print('Yielding dataset')
            yield dataset, coordinates, values

        except Exception as e:
            print(f'Unable to read line dataset {dataset}: {e}')
        finally:
            try:
                del netcdf_line_utils
            except:
                pass


def get_points_from_dict(dataset_point_dict):
    '''
    @param dataset_point_dict: {<dataset_path>: (<coordinates>, <values>),...}
    '''
    all_coordinates = []
    all_values = []
    for dataset in sorted(dataset_point_dict.keys()):
        coordinates, values = dataset_point_dict[dataset]
        all_coordinates += list(coordinates)
        all_values += list(values)

    print("Converting lists to arrays")
    all_values = np.array(all_values)
    all_coordinates = np.array(all_coordinates)
    assert all_values.shape[0] == all_coordinates.shape[0], 'Mismatched coordinate and value counts'
    print(f"A total of {all_values.shape[0]} points were read from {len(dataset_point_dict)} line datasets")

    return all_coordinates, all_values


def get_points_from_datasets(dataset_list,
                             wgs84_bounds,
                             variable_name,
                             coordinate_wkt,
                             min_points=None,
                             max_points=None):
    all_coordinates = []
    all_values = []
    for dataset, coordinates, values in read_datasets(dataset_list,
                                                      wgs84_bounds,
                                                      variable_name,
                                                      coordinate_wkt,
                                                      min_points,
                                                      max_points):
        all_coordinates += list(transform_coords(netcdf_line_utils.xycoords[spatial_selection_indices],
                                                 netcdf_line_utils.wkt,
                                                 coordinate_wkt))
        all_values += list(values)

    print("Converting lists to arrays")
    all_values = np.array(all_values)
    all_coordinates = np.array(all_coordinates)
    assert all_values.shape[0] == all_coordinates.shape[0], 'Mismatched coordinate and value counts'
    print(f"A total of {all_values.shape[0]} points were read from {len(dataset_list)} line datasets")

    return all_coordinates, all_values


def grid_points(coordinates,
                coordinate_wkt,
                values,
                grid_wkt,
                grid_bounds,
                grid_resolution,
                resampling_method='linear',
                point_step=1):
    '''
    Return interpolated grid from supplied coordinates and points
    '''

    # Determine spatial grid bounds rounded out to nearest GRID_RESOLUTION multiple
    pixel_centre_bounds = (round(math.floor(grid_bounds[0] / grid_resolution) * grid_resolution, 6),
                           round(math.floor(grid_bounds[1] / grid_resolution) * grid_resolution, 6),
                           round(math.floor(grid_bounds[2] / grid_resolution - 1.0) * grid_resolution + grid_resolution,
                                 6),
                           round(math.floor(grid_bounds[3] / grid_resolution - 1.0) * grid_resolution + grid_resolution,
                                 6)
                           )

    print("Reprojecting coordinates")
    grid_coordinates = np.array(transform_coords(coordinates, coordinate_wkt, grid_wkt))

    print("Computing spatial mask")
    spatial_subset_mask = np.logical_and(np.logical_and((grid_bounds[0] <= grid_coordinates[:, 0]),
                                                        (grid_coordinates[:, 0] <= grid_bounds[2])),
                                         np.logical_and((grid_bounds[1] <= grid_coordinates[:, 1]),
                                                        (grid_coordinates[:, 1] <= grid_bounds[3]))
                                         )
    # Create grids of Y and X values. Note YX ordering and inverted Y for image
    # Note GRID_RESOLUTION/2.0 fudge to avoid truncation due to rounding error
    print("Generating grid coordinates")
    grid_y, grid_x = np.mgrid[pixel_centre_bounds[3]:pixel_centre_bounds[1] - grid_resolution / 2.0:-grid_resolution,
                     pixel_centre_bounds[0]:pixel_centre_bounds[2] + grid_resolution / 2.0:grid_resolution]

    # Skip points to reduce memory requirements
    print("Generating spatial subset mask")
    point_subset_mask = np.zeros(shape=values.shape, dtype=bool)
    point_subset_mask[0:-1:point_step] = True
    point_subset_mask = np.logical_and(spatial_subset_mask, point_subset_mask)
    assert point_subset_mask.any(), f'No points found within grid bounds {grid_bounds}'

    grid_coordinates = grid_coordinates[point_subset_mask]

    # Interpolate required values to the grid - Note yx ordering for image
    print(f"Interpolating {grid_coordinates.shape[0]} points")
    grid_array = griddata(grid_coordinates[:, ::-1],
                          values[point_subset_mask],
                          (grid_y, grid_x),
                          method=resampling_method
                          )

    print("Interpolation complete")
    #  crs:GeoTransform = "109.1002342895272 0.00833333 0 -9.354948067227777 0 -0.00833333 "
    geotransform = [pixel_centre_bounds[0] - grid_resolution / 2.0,
                    grid_resolution,
                    0,
                    pixel_centre_bounds[3] + grid_resolution / 2.0,
                    0,
                    -grid_resolution
                    ]

    return grid_array, grid_wkt, geotransform
