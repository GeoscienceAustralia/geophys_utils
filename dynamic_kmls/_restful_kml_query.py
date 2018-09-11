'''
Created on 7 Sep. 2018

@author: Andrew Turner
'''
import yaml
import os
from flask_restful import Resource
from flask import request, make_response
import simplekml
import time
from shapely.geometry import Polygon
from shapely import wkt
from dynamic_kmls import netcdf2kml,  DEBUG, DATABASE_ENGINE
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache
import logging
#from pprint import pformat

# Define maximum bounding box width for point display. Uses survey convex-hull polygons for anything larger.
MAX_BOX_WIDTH_FOR_POINTS = 1.5

logger = logging.getLogger(__name__)
if DEBUG:
    logger.setLevel(logging.DEBUG) # Initial logging level for this module
else:
    logger.setLevel(logging.INFO) # Initial logging level for this module


settings = yaml.safe_load(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                            'netcdf2kml_settings.yml')))
#print('settings: {}'.format(pformat(settings)))
#print('settings: {}'.format(yaml.safe_dump(settings)))

class RestfulKMLQuery(Resource):
    '''
    RestfulKMLQuery Resource subclass for RESTful API
    '''
    def __init__(self):
        '''
        RestfulKMLQuery Constructor
        '''
        super(RestfulKMLQuery, self).__init__()
        
        self.sdmc = get_dataset_metadata_cache(db_engine=DATABASE_ENGINE, debug=DEBUG)
           
    
    def get(self, dataset_type):
        '''
        get method for RESTful API
        '''
        logger.debug('dataset_type: {}'.format(dataset_type))
        
        get_kml_functions = {'point': RestfulKMLQuery.build_point_dataset_kml,
                             'line': RestfulKMLQuery.build_line_dataset_kml,
                             'grid': RestfulKMLQuery.build_grid_dataset_kml
                             }
        
        dataset_settings = settings.get(dataset_type)
        
        #TODO: Handle this more gracefully
        assert dataset_settings, 'Invalid dataset type "{}"'.format(dataset_type)
        
        # Determine which KML generation function to call based on dataset_settings['format']
        get_kml_function = get_kml_functions.get(dataset_settings['format'])
        
        bbox = request.args['BBOX'] 
    
        xml = get_kml_function(self, bbox, dataset_type, dataset_settings)
        
        response = make_response(xml)
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
        

    def build_grid_dataset_kml(self, bbox, dataset_type, dataset_settings):
        '''
        Function to build KML for grid datasets
        @param bbox: Bounding box query parameter, e.g: BBOX=133.8666248233259,-16.80720537521252,135.0274640184073,-16.1150287021518
        @param dataset_type: dataset type string - must be a valid key in settings: e.g. 'aem' or 'ground_gravity'
        @param dataset_settings: settings for dataset_type as read from netcdf2kml_settings.yml settings file
        '''
        #TODO: Implement this
        pass
    
    
    def build_line_dataset_kml(self, bbox, dataset_type, dataset_settings):
        '''
        Function to build KML for line datasets
        @param bbox: Bounding box query parameter, e.g: BBOX=133.8666248233259,-16.80720537521252,135.0274640184073,-16.1150287021518
        @param dataset_type: dataset type string - must be a valid key in settings: e.g. 'aem' or 'ground_gravity'
        @param dataset_settings: settings for dataset_type as read from netcdf2kml_settings.yml settings file
        '''
        t0 = time.time()  # retrieve coordinates from query

        bbox_list = bbox.split(',')
        west = float(bbox_list[0])
        south = float(bbox_list[1])
        east = float(bbox_list[2])
        north = float(bbox_list[3])
    
        bbox_polygon = Polygon(((west, south),
                                (east, south),
                                (east, north),
                                (west, north),
                                (west, south)
                                ))
    
        t1 = time.time()
        logger.debug("Retrieve bbox values from get request...")
        logger.debug("Time: " + str(t1 - t0))
        
        logger.debug('Getting {} lines'.format(dataset_type))
    
        # Get the dataset_metadata_dict surveys from the database that are within the bbox
        
        dataset_metadata_dict_list = self.sdmc.search_dataset_distributions(
            keyword_list=dataset_settings['keyword_list'],
            protocol=dataset_settings['protocol'],
            ll_ur_coords=[[west, south], [east, north]]
        )
    
        logger.debug([[west, south], [east, north]])
        t2 = time.time()
        logger.debug("Retrieve dataset_metadata_dict strings from database...")
        logger.debug("Time: " + str(t2 - t1))
    
        kml = simplekml.Kml()
        dataset_type_folder = kml.newfolder(name=dataset_settings['netcdf_file_folder_name'])
    
        t_polygon_1 = time.time()
    
        if len(dataset_metadata_dict_list) > 0:
    
            for dataset_metadata_dict in dataset_metadata_dict_list:
                logger.debug("dataset_metadata_dict: {}".format(dataset_metadata_dict))
                
                netcdf_file_folder = dataset_type_folder.newfolder(name=dataset_metadata_dict['dataset_title'])

                netcdf_path = self.modify_nc_path(dataset_settings['netcdf_path_prefix'], str(dataset_metadata_dict['distribution_url']))
                
                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(netcdf_path, dataset_settings, dataset_metadata_dict)
                t_polygon_2 = time.time()
                logger.debug("set style and create netcdf2kmlconverter instance from dataset_metadata_dict for polygon ...")
                logger.debug("Time: " + str(t_polygon_2 - t_polygon_1))

                try:
                    survey_polygon = wkt.loads(dataset_metadata_dict.get('convex_hull_polygon'))
                except Exception as e:
                    logger.debug('Invalid geometry for polygon {}: {}'.format(dataset_metadata_dict.get('convex_hull_polygon'), e))
                    continue  # Skip this polygon

                if survey_polygon.intersects(bbox_polygon):
                # if survey_polygon.within(bbox_polygon):
                # if not survey_polygon.contains(bbox_polygon):
                # if survey_polygon.centroid.within(bbox_polygon):
                # if not survey_polygon.contains(bbox_polygon) and survey_polygon.centroid.within(bbox_polygon):

                    netcdf2kml_obj.build_lines(netcdf_file_folder, bbox_list)


                #dataset_polygon_region = netcdf2kml_obj.build_region(-1, -1, 200, 800)
            return str(netcdf_file_folder)
    
        else:
            empty_folder = kml.newfolder(name="No {} in view".format(dataset_settings['netcdf_file_folder_name']))
            return str(empty_folder)

    #grav
    
    #@app.route('/ground_gravity/<bounding_box>', methods=['GET'])
    def build_point_dataset_kml(self, bbox, dataset_type, dataset_settings):
        '''
        Function to build KML for point datasets
        @param bbox: Bounding box query parameter, e.g: BBOX=133.8666248233259,-16.80720537521252,135.0274640184073,-16.1150287021518
        @param dataset_type: dataset type string - must be a valid key in settings: e.g. 'aem' or 'ground_gravity'
        @param dataset_settings: settings for dataset_type as read from netcdf2kml_settings.yml settings file
        '''
        t0 = time.time()  # retrieve coordinates from query
    
        bbox_list = bbox.split(',')
        west = float(bbox_list[0])
        south = float(bbox_list[1])
        east = float(bbox_list[2])
        north = float(bbox_list[3])
    
        bbox_polygon = Polygon(((west, south),
                                (east, south),
                                (east, north),
                                (west, north),
                                (west, south)
                                ))
    
        t1 = time.time()
        logger.debug("Retrieve bbox values from get request...")
        logger.debug("Time: " + str(t1 - t0))
    
        # Get the dataset_metadata_dict surveys from the database that are within the bbox
        dataset_metadata_dict_list = self.sdmc.search_dataset_distributions(
            keyword_list=dataset_settings['keyword_list'],
            protocol=dataset_settings['protocol'],
            ll_ur_coords=[[west, south], [east, north]]
        )
        logger.debug("tuple: " + str(dataset_metadata_dict_list))
    
        logger.debug([[west, south], [east, north]])
        t2 = time.time()
        logger.debug("Retrieve dataset_metadata_dict strings from database...")
        logger.debug("Time: " + str(t2 - t1))
    
        kml = simplekml.Kml()
        netcdf_file_folder = kml.newfolder(name=dataset_settings['netcdf_file_folder_name'])
    
        # ----------------------------------------------------------------------------------------------------------------
        # High zoom: show points rather than polygons.
        if east - west < MAX_BOX_WIDTH_FOR_POINTS:
            logger.debug('Getting {} points'.format(dataset_type))
            if len(dataset_metadata_dict_list) > 0:
                for dataset_metadata_dict in dataset_metadata_dict_list:                   
                    netcdf_path = self.modify_nc_path(dataset_settings['netcdf_path_prefix'], str(dataset_metadata_dict['distribution_url']))
                    
                    logger.debug("Building NETCDF: {} ".format(netcdf_path))
                    netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(netcdf_path, dataset_settings, dataset_metadata_dict)
                    t3 = time.time()
                    logger.debug("set style and create netcdf2kmlconverter instance...")
                    logger.debug("Time: " + str(t3 - t2))
    
                    # logger.debug("Number of points in file: " + str(netcdf2kml_obj.npu.point_count))
    
                    ta = time.time()
                    netcdf2kml_obj.build_points(netcdf_file_folder, bbox_list)
                    tb = time.time()
                    logger.debug("do the things time: " + str(tb - ta))
                    logger.debug("Build the point ...")
                    dataset_points_region = netcdf2kml_obj.build_region(100, -1, 200, 800)
                    netcdf_file_folder.region = dataset_points_region
                    netcdf2kml_obj.netcdf_dataset.close()  # file must be closed after use to avoid errors when accessed again.
                    del netcdf2kml_obj  # Delete netcdf2kml_obj to removenetcdf2kml_obj.npu cache file
                    t4 = time.time()
                return str(netcdf_file_folder)
    
            else:
                empty_folder = kml.newfolder(name="No {} in view".format(dataset_settings['netcdf_file_folder_name']))
                return empty_folder
    
        # ----------------------------------------------------------------------------------------------------------------
        # Low zoom: show polygons and not points.
        else:
            logger.debug('Getting {} polygons'.format(dataset_type))
            t_polygon_1 = time.time()
    
            if len(dataset_metadata_dict_list) > 0:
    
                for dataset_metadata_dict in dataset_metadata_dict_list:
                    netcdf_path = self.modify_nc_path(dataset_settings['netcdf_path_prefix'], str(dataset_metadata_dict['distribution_url']))
                    #logger.debug(netcdf_path)
                    netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(netcdf_path, dataset_settings, dataset_metadata_dict)
                    t_polygon_2 = time.time()
                    logger.debug("set style and create netcdf2kmlconverter instance from dataset_metadata_dict for polygon ...")
                    logger.debug("Time: " + str(t_polygon_2 - t_polygon_1))
    
                    try:
                        survey_polygon = wkt.loads(dataset_metadata_dict['convex_hull_polygon'])
                    except Exception as e:
                        # print(e)
                        continue  # Skip this polygon
    
                    if survey_polygon.intersects(bbox_polygon):
                    # if survey_polygon.within(bbox_polygon):
                    # if not survey_polygon.contains(bbox_polygon):
                    # if survey_polygon.centroid.within(bbox_polygon):
                    # if not survey_polygon.contains(bbox_polygon) and survey_polygon.centroid.within(bbox_polygon):

                        polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder)
                        #polygon_folder = netcdf2kml_obj.build_lines(netcdf_file_folder, bbox_list)

    
                    dataset_polygon_region = netcdf2kml_obj.build_region(-1, -1, 200, 800)
                # polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
    
                # else:  # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                #    logger.debug("not enough points")
    
                # neww = kml.save("test_polygon.kml")
                return str(netcdf_file_folder)
            else:
                empty_folder = kml.newfolder(name="No {} in view".format(dataset_settings['netcdf_file_folder_name']))
                return str(empty_folder)