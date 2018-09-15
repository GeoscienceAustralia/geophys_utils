'''
Created on 7 Sep. 2018

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
import os
from flask_restful import Resource
from flask import request, make_response
import simplekml
from shapely.geometry import Polygon
from dynamic_kmls.netcdf2kml import NetCDF2kmlConverter, settings
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache
import logging
#from pprint import pformat

logger = logging.getLogger(__name__)
if settings['global_settings']['debug']:
    logger.setLevel(logging.DEBUG) # Initial logging level for this module
else:
    logger.setLevel(logging.INFO) # Initial logging level for this module

class RestfulKMLQuery(Resource):
    '''
    RestfulKMLQuery Resource subclass for RESTful API
    '''
    def __init__(self):
        '''
        RestfulKMLQuery Constructor
        '''
        super(RestfulKMLQuery, self).__init__()
        
        self.sdmc = get_dataset_metadata_cache(db_engine=settings['global_settings']['database_engine'], 
                                               debug=settings['global_settings']['debug'])
           
    
    def get(self, dataset_type):
        '''
        get method for RESTful API
        '''
        logger.debug('dataset_type: {}'.format(dataset_type))
        
        dataset_settings = settings['dataset_settings'].get(dataset_type)
        
        #TODO: Handle this more gracefully
        assert dataset_settings, 'Invalid dataset type "{}"'.format(dataset_type)
        
        bbox_list = [float(ordinate) for ordinate in request.args['BBOX'].split(',')]
    
        #=======================================================================
        # # return polygon KML for low zoom
        # if ((bbox_list[2] - bbox_list[0]) >= dataset_settings['min_polygon_bbox_width']):
        #     kml = self.build_polygon_kml(bbox_list, dataset_type, settings)
        # else: # High zoom - plot detailed dataset
        #     kml = get_kml_function(self, bbox_list, dataset_type, settings)
        #=======================================================================
        
        dataset_metadata_dict_list = self.sdmc.search_dataset_distributions(
            keyword_list=dataset_settings['keyword_list'],
            protocol=dataset_settings['protocol'],
            ll_ur_coords=[[bbox_list[0], bbox_list[1]], [bbox_list[2], bbox_list[3]]]
        )
    
        kml = simplekml.Kml()
    
        if len(dataset_metadata_dict_list) > 0:
            dataset_type_folder = kml.newfolder(name=dataset_settings['dataset_type_name'])
    
            for dataset_metadata_dict in dataset_metadata_dict_list:
                logger.debug("dataset_metadata_dict: {}".format(dataset_metadata_dict))
                
                # dataset_folder = dataset_type_folder.newfolder(name=dataset_metadata_dict['dataset_title'])

                #TODO: Replace this hack for turning OPeNDAP endpoints into local pathnames
                netcdf_path = self.modify_nc_path(dataset_settings['netcdf_path_prefix'], 
                                                  str(dataset_metadata_dict['distribution_url']))
                
                netcdf2kml_obj = NetCDF2kmlConverter(netcdf_path, settings, dataset_type, dataset_metadata_dict)
                netcdf2kml_obj.build_kml(dataset_type_folder, bbox_list)

                del netcdf2kml_obj # Finished with this object - delete it
    
        else:
            dataset_type_folder = kml.newfolder(name="No {} in view".format(dataset_settings['dataset_type_name']))
        
        response = make_response(str(dataset_type_folder))
        response.headers['content-type'] = 'application/vnd.google-earth.kml+xml'
        return response

    
    def modify_nc_path(self, netcdf_path_prefix, opendap_endpoint):    
        '''
        Helper function to substitute netcdf_path_prefix in netcdf_path if defined
        '''
        if netcdf_path_prefix:
            return os.path.join(netcdf_path_prefix, os.path.basename(opendap_endpoint))
        else:
            return opendap_endpoint
        

    def rectangle_from_bbox_list(self, bbox_list):
        '''
        Helper function to return bounding box rectangle from [<xmin>, <ymin>, <xmax>, <ymax>] list
        '''
        return Polygon(((bbox_list[0], bbox_list[1]),
                        (bbox_list[2], bbox_list[1]),
                        (bbox_list[2], bbox_list[3]),
                        (bbox_list[0], bbox_list[3]),
                        (bbox_list[0], bbox_list[1])
                        ))
        
