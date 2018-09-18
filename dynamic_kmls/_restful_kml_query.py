'''
Created on 7 Sep. 2018

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
import os
from flask_restful import Resource
from flask import request, make_response
from shapely.geometry import Polygon
from dynamic_kmls import settings
from dynamic_kmls.netcdf2kml import NetCDF2kmlConverter
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache
import logging
#from pprint import pformat

logger = logging.getLogger(__name__)
    
if settings['global_settings']['debug']:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)
logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))

class RestfulKMLQuery(Resource):
    '''
    RestfulKMLQuery Resource subclass for RESTful API
    '''
    CONTENT_TYPE = 'application/vnd.google-earth.kml+xml'
    
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
        #  logger.debug("dataset_metadata_dict_list: {}".format(dataset_metadata_dict_list))
            
        # Insert modified netCDF file path into each dict in list
        for dataset_metadata_dict in dataset_metadata_dict_list:
            #TODO: Replace this hack for turning OPeNDAP endpoints into local pathnames
            dataset_metadata_dict['netcdf_path'] = self.modify_nc_path(dataset_settings['netcdf_path_prefix'], 
                                                                       str(dataset_metadata_dict['distribution_url']))
            
        netcdf2kml_obj = NetCDF2kmlConverter(settings, 
                                             dataset_type, 
                                             request_host = request.host,
                                             debug=settings['global_settings']['debug'])        
        netcdf2kml_obj.build_bbox_kml(dataset_metadata_dict_list, bbox_list)
   
        response = make_response(netcdf2kml_obj.kml_string)
        response.headers['content-type'] = RestfulKMLQuery.CONTENT_TYPE
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
        

