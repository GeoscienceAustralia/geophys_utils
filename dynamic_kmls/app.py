from flask import Flask
from flask_restful import Api
from flask import request
import simplekml
import netCDF4
import time
import os
import re
from shapely.geometry import Polygon
from shapely import wkt
from geophys_utils import NetCDFPointUtils
from geophys_utils.netcdf_converter import netcdf2kml
from geophys_utils.dataset_metadata_cache import SQLiteDatasetMetadataCache
import logging

# Define maximum bounding box width for point display. Uses survey convex-hull polygons for anything larger.
MAX_BOX_WIDTH_FOR_POINTS = 1.5

# Set the following to None or empty string to use OPeNDAP endpoints
LOCAL_FILE_LOCATION = None
#LOCAL_FILE_LOCATION = 'D:\Temp\gravity point_datasets'

app = Flask(__name__)
api = Api(app)

# Setup logging handlers if required
logger = logging.getLogger(__name__)  # Get logger
logger.setLevel(logging.DEBUG)  # Initial logging level for this module


@app.route('/<bounding_box>',methods=['GET'])
def do_everything(bounding_box):

    t0 = time.time()  # retrieve coordinates from query

    query_string = request.args
    bbox = query_string['BBOX']
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
    logger.debug("Time: " + str(t1-t0))

    # Get the point_data_tuple surveys from the database that are within the bbox
    sdmc = SQLiteDatasetMetadataCache(debug=False)
    point_data_tuple_list = sdmc.search_dataset_distributions(
        keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
        protocol='opendap',
        ll_ur_coords=[[west, south], [east, north]]
        )

    logger.debug([[west, south], [east, north]])
    t2 = time.time()
    logger.debug("Retrieve point_data_tuple strings from database...")
    logger.debug("Time: " + str(t2-t1))

    kml = simplekml.Kml()

    # ----------------------------------------------------------------------------------------------------------------
    # High zoom: show points rather than polygons.
    if east - west < MAX_BOX_WIDTH_FOR_POINTS:

        if len(point_data_tuple_list) > 0:
            # set point style
            point_style = simplekml.Style()
            point_style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
            point_style.iconstyle.scale = 0.7
            point_style.labelstyle.scale = 0  # removes the label

            netcdf_file_folder = kml.newfolder(name="Ground Gravity Survey Observations")

            for point_data_tuple in point_data_tuple_list:
                logger.debug("Building NETCDF: " + str(point_data_tuple[2]))
                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(point_data_tuple)
                t3 = time.time()
                logger.debug("set style and create netcdf2kmlconverter instance of point_data_tuple file ...")
                logger.debug("Time: " + str(t3 - t2))

                #logger.debug("Number of points in file: " + str(netcdf2kml_obj.npu.point_count))

                nc_path = netcdf2kml_obj.netcdf_path
                if LOCAL_FILE_LOCATION:
                    nc_path = os.path.join(LOCAL_FILE_LOCATION,
                                           os.path.basename(nc_path)
                                           )
                    
                netcdf2kml_obj.netcdf_dataset = netCDF4.Dataset(nc_path)

                netcdf2kml_obj.npu = NetCDFPointUtils(netcdf2kml_obj.netcdf_dataset)

                if netcdf2kml_obj.npu.point_count > 0:
                    ta = time.time()
                    netcdf2kml_obj.build_points(netcdf_file_folder, bbox_list, point_style)
                    tb = time.time()
                    logger.debug("do the things time: " + str(tb-ta))
                    logger.debug("Build the point ...")
                dataset_points_region = netcdf2kml_obj.build_region(100, -1, 200, 800)
                netcdf_file_folder.region = dataset_points_region
                netcdf2kml_obj.netcdf_dataset.close() # file must be closed after use to avoid errors when accessed again.
                del netcdf2kml_obj # Delete netcdf2kml_obj to removenetcdf2kml_obj.npu cache file
                t4 = time.time()

            return str(netcdf_file_folder)

        else:
            logger.debug("No surveys in view")

    # ----------------------------------------------------------------------------------------------------------------
    # Low zoom: show polygons and not points.
    else:
        t_polygon_1 = time.time()

        # set polygon style
        polygon_style = simplekml.Style()
        polygon_style.polystyle.color = 'B30000ff'  # Transparent red
        #polygon_style.polystyle.color = 'ff4545'
        polygon_style.polystyle.outline = 1

        polygon_style_background = simplekml.Style()
        polygon_style_background.polystyle.color = '7FFFFFFF'  # Transparent white
        polygon_style_background.polystyle.outline = 1

        if len(point_data_tuple_list) > 0:
            netcdf_file_folder = kml.newfolder(name="Ground Gravity Survey Extents")
            for point_data_tuple in point_data_tuple_list:
                logger.debug("NETCDF: " + str(point_data_tuple))
            
                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(point_data_tuple)
                t_polygon_2 = time.time()
                logger.debug("set style and create netcdf2kmlconverter instance of point_data_tuple file for polygon ...")
                logger.debug("Time: " + str(t_polygon_2 - t_polygon_1))
            
                try:
                    survey_polygon = wkt.loads(point_data_tuple[3])
                except Exception as e:
                    #print(e)
                    continue # Skip this polygon
                    
                if survey_polygon.within(bbox_polygon):
                #if not survey_polygon.contains(bbox_polygon):
                #if survey_polygon.centroid.within(bbox_polygon):
                #if not survey_polygon.contains(bbox_polygon) and survey_polygon.centroid.within(bbox_polygon):
                    
                    polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder, polygon_style)
                else:
                    polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder, polygon_style, False)
            
                dataset_polygon_region = netcdf2kml_obj.build_region(-1, -1, 200, 800)
                polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
            
                #else:  # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                #    logger.debug("not enough points")
        
            # neww = kml.save("test_polygon.kml")
            return str(netcdf_file_folder)

        else:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)

app.run(debug=True)

