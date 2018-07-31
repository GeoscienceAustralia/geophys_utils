from flask import Flask
from flask_restful import Api, Resource, reqparse
from flask import request
import simplekml
from geophys_utils import NetCDFPointUtils
from geophys_utils import points2convex_hull
from geophys_utils import _polygon_utils
import netCDF4
import time
from geophys_utils import NetCDFPointUtils
from geophys_utils.netcdf_converter import netcdf2kml
import argparse
from geophys_utils.dataset_metadata_cache import SQLiteDatasetMetadataCache

app = Flask(__name__)
api = Api(app)

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



    t1 = time.time()
    print("Retrieve bbox values from get request...")
    print("Time: " + str(t1-t0))

    # Get the netcdf surveys from the database that are within the bbox
    sdmc = SQLiteDatasetMetadataCache(debug=True)
    endpoint_list = sdmc.search_dataset_distributions(
        keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
        protocol='opendap',
        ll_ur_coords=[[west, south], [east, north]]
        )
    print([[west, south], [east, north]])
    t2 = time.time()
    print("Retrieve netcdf strings from database...")
    print("Time: " + str(t2-t1))
    #print("ENDPOINT LIST" + str(endpoint_list))


    kml = simplekml.Kml()
    # ----------------------------------------------------------------------------------------------------------------
    # low zoom, show points only.
    if east - west < 1:

        if len(endpoint_list) > 0:

            # set point style
            point_style = simplekml.Style()
            point_style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
            point_style.iconstyle.scale = 0.7
            point_style.labelstyle.scale = 0  # removes the label

            netcdf_file_folder = kml.newfolder()

            for netcdf in endpoint_list:

                print("Building NETCDF: " + str(netcdf[2]))

                netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)

                t3 = time.time()
                print("set style and create netcdf2kmlconverter instance of netcdf file ...")
                print("Time: " + str(t3 - t2))

                #print("Number of points in file: " + str(netcdf2kml_obj.npu.point_count))


                dataset_points_region = netcdf2kml_obj.build_region(100, -1, 200, 800)
                netcdf_file_folder.region = dataset_points_region
                print("HERE")
                print(netcdf2kml_obj.survey_title)
                ta = time.time()
                new_survey_folder = netcdf_file_folder.newfolder(name=str(netcdf2kml_obj.survey_title) + " " +
                                                                          str(netcdf2kml_obj.survey_id))
                print(new_survey_folder)
                print(bbox_list)
                print(point_style)
                new_survey_folder = netcdf2kml_obj.build_points(new_survey_folder, bbox_list, point_style)

                tb = time.time()
                print("do the things time: " + str(tb-ta))

                #print(netcdf_file_folder)
                print("Build the point ...")
                #print("Time: " + str(t4 - t3))
                #print("ENDPOINT LIST" + str(endpoint_list))

            t4 = time.time()

            return str(netcdf_file_folder)

        else:
            print("No surveys in view")

    # ----------------------------------------------------------------------------------------------------------------
    # Zoomed in, show polygons only.
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

        if len(endpoint_list) > 0:
             netcdf_file_folder = kml.newfolder()
             for netcdf in endpoint_list:
                 print("NETCDF: " + str(netcdf))



                 netcdf2kml_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
                 t_polygon_2 = time.time()
                 print("set style and create netcdf2kmlconverter instance of netcdf file for polygon ...")
                 print("Time: " + str(t_polygon_2 - t_polygon_1))


                 #if netcdf2kml_obj.npu.point_count > 2:
                 #if netcdf[4]: # if gridflag is 1 (true)
                 if netcdf[8]:
                    polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder, polygon_style)
                 else:
                     polygon_folder = netcdf2kml_obj.build_polygon(netcdf_file_folder, polygon_style, False)

                 dataset_polygon_region = netcdf2kml_obj.build_region(-1, -1, 200, 800)
                 polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder

                 #else:  # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                  #  print("not enough points")

            # neww = kml.save("test_polygon.kml")
             return str(netcdf_file_folder)

        else:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)

app.run(debug=True)

