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
    # print('query done!')
    # print(time.time())
    bbox = query_string['BBOX']
    bbox_list = bbox.split(',')
    west = float(bbox_list[0]) + 4
    south = float(bbox_list[1]) - 4
    east = float(bbox_list[2]) - 4
    north = float(bbox_list[3]) + 4


    if east - west < 5:
        print("GET POINTS")
        sdmc = SQLiteDatasetMetadataCache(debug=True)

        endpoint_list = sdmc.search_dataset_distributions(keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
                                                     protocol='opendap',
                                                     #ll_ur_coords=[[-179.9, -90.0], [180.0, 90.0]]
                                                          # if search with no paramters reutrns everything.
                                                          #ll lower left, upper right
                                                     ll_ur_coords=[[west, south], [east, north]]
                                                     )
        print("ENDPOINT LIST" + str(endpoint_list))
        kml = simplekml.Kml()

        if len(endpoint_list) > 0:
            netcdf_file_folder = kml.newfolder()
            all_points_folder = kml.newfolder()
            for netcdf in endpoint_list:


                print("NETCDF: " + str(netcdf))

                converter_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
                print(converter_obj.npu.point_count)
                if converter_obj.npu.point_count > 2:

                    polygon_folder, polygon = converter_obj.build_polygon(netcdf_file_folder)

                    dataset_polygon_region = converter_obj.build_region(converter_obj.MAX_LOD_PIXELS, -1,
                                                               converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    dataset_points_region = converter_obj.build_region(converter_obj.MIN_LOD_PIXELS - 200, converter_obj.MAX_LOD_PIXELS,
                                                                       converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
                    survey_points_folder = kml.newfolder()
                    #bbox_list = [west, south, east, north]
                    #bbox_list = [west, south, east, north]
                    #bbox_list = [west, south, east, north]
                    #bbox_list = [south, east, north, west]
                    survey_points_folder = converter_obj.do_the_things(survey_points_folder, kml, bbox_list)
                    print(survey_points_folder)
                    #all_points_folder = all_points_folder.survey_points_folder
                    all_points_folder.region = dataset_points_region  # insert built point region into point folder
                    #points_folder = converter_obj.do_the_things(netcdf, points_folder, query_string)

                    #all_the_points = all_the_points + str(points_folder)

                elif converter_obj.npu.point_count == 0:
                    #empty_folder = kml.newfolder(name="no points in view")
                    #return str(empty_folder)
                    print('nada')

                # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                else:
                    print("not enough points")
                    #empty_folder = kml.newfolder(name="no points in view")
                    #return str(empty_folder)

            new = kml.save("testt.kml")
            print("saved: " + str(new))
            return str(netcdf_file_folder)
                #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))

        # then add the network link using this args.server
        elif len(endpoint_list) == 0:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)

        # for surveys with 1 or 2 points. Can't make a polygon.
        else:
            print("not enough points")
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)
            # single
            #print("single survey")


    else:  # polygons only

        print("anything")

        sdmc = SQLiteDatasetMetadataCache(debug=True)

        endpoint_list = sdmc.search_dataset_distributions(keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
                                                     protocol='opendap',
                                                     #ll_ur_coords=[[-179.9, -90.0], [180.0, 90.0]]
                                                          # if search with no paramters reutrns everything.
                                                          #ll lower left, upper right
                                                     ll_ur_coords=[[west, south], [east, north]]
                                                     )
        print("ENDPOINT LIST" + str(endpoint_list))
        # query_string = request.args
        # parser = argparse.ArgumentParser()
        #
        # # parser.add_argument("-s", "--server",
        # #                     help="The server to receive the get request from google earth and dynamically "
        # #                          "build kml points within the bbox. If this parameter is empty, a static "
        # #                          "kml will be generated", type=str, required=False)
        # parser.add_argument("-n", "--netcdf_path_list", help="Add one or more paths to netcdf files to be converted into a"
        #                                                      "single kml file.", type=str, nargs='+')
        #
        # args = parser.parse_args()

        kml = simplekml.Kml()


        # if len(args.netcdf_path_list) > 1:
        #     print("multiple surveys")
        #     # multiples
           # list_of_point_folders= []

        if len(endpoint_list) > 0:
            all_the_points = ""
            netcdf_file_folder = kml.newfolder()
            for netcdf in endpoint_list:
                print("NETCDF: " + str(netcdf))
                #points_folder = kml.newfolder()

                converter_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
                print(converter_obj.npu.point_count)
                if converter_obj.npu.point_count > 2:

                    polygon_folder, polygon = converter_obj.build_polygon(netcdf_file_folder)

                    dataset_polygon_region = converter_obj.build_region(converter_obj.MAX_LOD_PIXELS, -1,
                                                               converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    #dataset_points_region = converter_obj.build_region(converter_obj.MIN_LOD_PIXELS - 200, converter_obj.MAX_LOD_PIXELS,
                     #                                                  converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder

                    #points_folder.region = dataset_points_region  # insert built point region into point folder
                    #points_folder = converter_obj.do_the_things(netcdf, points_folder, query_string)
                    #points_folder = do_the_things(netcdf, points_folder)
                    #all_the_points = all_the_points + str(points_folder)

                elif converter_obj.npu.point_count == 0:
                    #empty_folder = kml.newfolder(name="no points in view")
                    #return str(empty_folder)
                    print('nada')

                # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                else:
                    print("not enough points")
                    #empty_folder = kml.newfolder(name="no points in view")
                    #return str(empty_folder)

            new = kml.save("testt.kml")
            print("saved: " + str(new))
            return str(netcdf_file_folder)
                #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))

        # then add the network link using this args.server
        elif len(endpoint_list) == 0:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)

        # for surveys with 1 or 2 points. Can't make a polygon.
        else:
            print("not enough points")
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)
            # single
            #print("single survey")

            #folder = do_the_things(endpoint_list[0], points_folder)
            # converter_object = netcdf2kml.NetCDF2kmlConverter(args.netcdf_path_list[0])
            # converter_object.build_dynamic_kml()
            # converter_object.kml.save(converter_object.survey_title + " dynamic points.kml")



app.run(debug=True)

