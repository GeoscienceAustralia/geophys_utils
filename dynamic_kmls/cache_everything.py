'''
Created on 3 Oct. 2018

@author: Alex Ip - Geoscience Australia

Hacky utility to cache EVERYTHING from OPeNDAP & WMS endpoints. 
WARNING: This will take some time to run to completion and consume a considerable amount of bandwidth. Do NOT try this at home.
'''
import os
import tempfile
from dynamic_kmls import settings
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils
from dynamic_kmls import cache_image_file

bbox_list = [-180.0, -90.0, 180.0, 90.0] # Get everything regardless of spatial position

netcdf_util_subclass={'point': NetCDFPointUtils, 'line': NetCDFLineUtils} # Subclasses keyed by dataset_format

def main():
    dataset_metadata_cache = get_dataset_metadata_cache(db_engine=settings['global_settings']['database_engine'], 
                                                        debug=settings['global_settings']['debug'])
    
    for dataset_type, dataset_settings in settings['dataset_settings'].items():
        
        dataset_format = dataset_settings['dataset_format']
        
        if dataset_format not in ['point', 'line', 'grid']:
            continue
        
        cache_dir = os.path.join((settings['global_settings'].get('cache_root_dir') or 
                          tempfile.gettempdir()),
                          'kml_server_cache',
                          dataset_type
                          )
        os.makedirs(cache_dir, exist_ok=True)
        
        dataset_metadata_dict_list = dataset_metadata_cache.search_dataset_distributions(
            keyword_list=dataset_settings['keyword_list'],
            protocol=dataset_settings['protocol'],
            ll_ur_coords=[[bbox_list[0], bbox_list[1]], [bbox_list[2], bbox_list[3]]]
        )
        print('{} {} datasets found'.format(len(dataset_metadata_dict_list), dataset_type))
        
        
        for dataset_metadata_dict in dataset_metadata_dict_list:
            distribution_url = dataset_metadata_dict['distribution_url']
            
            if dataset_format == 'grid':
                wms_url = distribution_url.replace('/dodsC/', '/wms/') #TODO: Replace this hack
                
                image_basename = os.path.splitext(os.path.basename(distribution_url))[0] + '.png'
                
                wms_pixel_size = dataset_settings.get('wms_pixel_size') or settings['default_dataset_settings'].get('wms_pixel_size')
                assert wms_pixel_size, 'Unable to determine wms_pixel_size'
                
                # Retrieve image for entire dataset
                north = dataset_metadata_dict['latitude_max']
                south = dataset_metadata_dict['latitude_min']
                east = dataset_metadata_dict['longitude_max']
                west = dataset_metadata_dict['longitude_min']

                wms_url = wms_url + "?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&BBOX={0},{1},{2},{3}&CRS=EPSG:4326&WIDTH={4}&HEIGHT={5}&LAYERS={6}&STYLES=&FORMAT=image/png" \
                      "&DPI=120&MAP_RESOLUTION=120&FORMAT_OPTIONS=dpi:120&TRANSPARENT=TRUE" \
                      "&COLORSCALERANGE={7}%2C{8}&NUMCOLORBANDS=127".format(south, 
                                                                             west, 
                                                                             north, 
                                                                             east, 
                                                                             int((east - west) / wms_pixel_size),
                                                                             int((north - south) / wms_pixel_size),
                                                                             dataset_settings['wms_layer_name'],
                                                                             dataset_settings['wms_color_range'][0],
                                                                             dataset_settings['wms_color_range'][1],
                                                                             )
                      
                print('\tCaching image {} from {}'.format(image_basename, wms_url))

                cached_image_url_path = cache_image_file(dataset_type, image_basename, wms_url)
                
                print('\t\tImage URL: {}'.format(cached_image_url_path))
                
                continue
            
            # Points and lines handled below
            
            try:
                print('\tOpening {}'.format(distribution_url))
                
                netcdf_util = netcdf_util_subclass[dataset_format](distribution_url, 
                     enable_disk_cache=True,
                     enable_memory_cache=True,
                     cache_dir=cache_dir,
                     debug=settings['global_settings']['debug'])
                
                print('\t\tCached {} points'.format(len(netcdf_util.xycoords))) # Cause xycoords to be cached
                
                if dataset_type == 'line':
                    print('\t\tCached {} lines'.format(len(netcdf_util.line))) # Cause line & line_index to be cached
                    
                netcdf_util.netcdf_dataset.close()
            except BaseException as e:
                print('\t\tUnable to cache data for {}: {}'.format(distribution_url, e))
    

if __name__ == '__main__':
    main()