from flask import Flask
from flask_restful import Api
from flask import request
import simplekml
import time
from shapely.geometry import Polygon
from shapely import wkt
from geophys_utils.netcdf_converter import netcdf2kml
from geophys_utils.dataset_metadata_cache import get_dataset_metadata_cache
import logging
import argparse
import yaml
import os

DATABASE_ENGINE = 'SQLite'
# DATABASE_ENGINE = 'Postgres'

# Define maximum bounding box width for point display. Uses survey convex-hull polygons for anything larger.
MAX_BOX_WIDTH_FOR_POINTS = 1.5

app = Flask(__name__)
api = Api(app)

# Setup logging handlers if required
logger = logging.getLogger(__name__)  # Get logger
logger.setLevel(logging.DEBUG)  # Initial logging level for this module

# while i < len(sys.argv):
#     print(sys.argv[i])
#     netcdf_path_list.append(sys.argv[i])
#     i += 1
# print(netcdf_path_list)

parser = argparse.ArgumentParser()
#
parser.add_argument("-d", "--dataset_settings", help="Point the flask server to the correct dataset settings in the "
                                                     "yaml file. Options include: ground_gravity, aem.", type=str,
                    required=False)

args = parser.parse_args()

settings = yaml.safe_load(open(os.path.dirname(__file__) + '/netcdf2kml_settings.yml'))
yaml_settings = settings[args.dataset_settings]


@app.route('/<bounding_box>', methods=['GET'])
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
    logger.debug("Time: " + str(t1 - t0))

    # Get the point_data_tuple surveys from the database that are within the bbox
    sdmc = get_dataset_metadata_cache(db_engine=DATABASE_ENGINE, debug=False)
    point_data_tuple_list = sdmc.search_dataset_distributions(
        keyword_list=yaml_settings['keyword_list'],
        protocol=yaml_settings['protocol'],
        ll_ur_coords=[[west, south], [east, north]]
    )

    logger.debug([[west, south], [east, north]])
    t2 = time.time()
    logger.debug("Retrieve point_data_tuple strings from database...")
    logger.debug("Time: " + str(t2 - t1))

    kml = simplekml.Kml()
    netcdf_file_folder = kml.newfolder(name=yaml_settings['netcdf_file_folder_name'])

    # ----------------------------------------------------------------------------------------------------------------
    # High zoom: show points rather than polygons.
    if east - west < MAX_BOX_WIDTH_FOR_POINTS:

        if len(point_data_tuple_list) > 0:
            for point_data_tuple in point_data_tuple_list:
                logger.debug("Building NETCDF: " + str(point_data_tuple[2]))
                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(yaml_settings, point_data_tuple)
                t3 = time.time()
                logger.debug("set style and create netcdf2kmlconverter instance of point_data_tuple file ...")
                logger.debug("Time: " + str(t3 - t2))

                # logger.debug("Number of points in file: " + str(netcdf2kml_obj.npu.point_count))

                if netcdf2kml_obj.npu.point_count > 0:
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
            logger.debug("No surveys in view")

    # ----------------------------------------------------------------------------------------------------------------
    # Low zoom: show polygons and not points.
    else:
        t_polygon_1 = time.time()

        if len(point_data_tuple_list) > 0:

            for point_data_tuple in point_data_tuple_list:
                logger.debug("point_data_tuple: " + str(point_data_tuple))
                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(yaml_settings, point_data_tuple)
                t_polygon_2 = time.time()
                logger.debug("set style and create netcdf2kmlconverter instance from point_data_tuple for polygon ...")
                logger.debug("Time: " + str(t_polygon_2 - t_polygon_1))

                try:
                    survey_polygon = wkt.loads(point_data_tuple[3])
                except Exception as e:
                    # print(e)
                    continue  # Skip this polygon

                if survey_polygon.intersects(bbox_polygon):
                    # if survey_polygon.within(bbox_polygon):
                    # if not survey_polygon.contains(bbox_polygon):
                    # if survey_polygon.centroid.within(bbox_polygon):
                    # if not survey_polygon.contains(bbox_polygon) and survey_polygon.centroid.within(bbox_polygon):

                    # polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder)
                    polygon_folder = netcdf2kml_obj.build_lines(netcdf_file_folder, bbox_list)

                else:
                    polygon_folder = netcdf2kml_obj.build_lines(netcdf_file_folder, False)

                dataset_polygon_region = netcdf2kml_obj.build_region(-1, -1, 200, 800)
            # polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder

            # else:  # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
            #    logger.debug("not enough points")

            # neww = kml.save("test_polygon.kml")
            return str(netcdf_file_folder)
        else:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)


app.run(debug=True)
