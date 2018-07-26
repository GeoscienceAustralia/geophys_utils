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
    west = float(bbox_list[0])
    south = float(bbox_list[1])
    east = float(bbox_list[2])
    north = float(bbox_list[3])

    print("GET POINTS")
    t1 = time.time()
    sdmc = SQLiteDatasetMetadataCache(debug=True)

    endpoint_list = sdmc.search_dataset_distributions(
        keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
        protocol='opendap',
        # ll_ur_coords=[[-179.9, -90.0], [180.0, 90.0]]
        # if search with no paramters reutrns everything.
        # ll lower left, upper right
        ll_ur_coords=[[west, south], [east, north]]
        )
    t2 = time.time()
    print("ENDPOINT LIST" + str(endpoint_list))
    kml = simplekml.Kml()

    # low zoom, show polygons only. High zoom show points
    if east - west < 1:

        if len(endpoint_list) > 0:
            netcdf_file_folder = kml.newfolder()
            for netcdf in endpoint_list:


                print("NETCDF: " + str(netcdf))

                converter_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
                print(converter_obj.npu.point_count)
                if converter_obj.npu.point_count > 2:

#                    polygon_folder, polygon = converter_obj.build_polygon(netcdf_file_folder)

                    #dataset_polygon_region = converter_obj.build_region(converter_obj.MAX_LOD_PIXELS, -1,
                                                           #    converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    dataset_points_region = converter_obj.build_region(100, -1,
                                                                       converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)

                    netcdf_file_folder.region = dataset_points_region


                    ta = time.time()
                    converter_obj.do_the_things(netcdf_file_folder, kml, bbox_list)


                    tb = time.time()
                    print("do the things time: " + str(tb-ta))



                    print(netcdf_file_folder)

                    #all_points_folder = all_points_folder.survey_points_folder
                    # insert built point region into point folder
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

            new = kml.save("testtt.kml")
            print("saved: " + str(new))
            t3 = time.time()
            print("query time: " + str(t1 - t0))
            print("database search time: " + str(t2 - t1))
            print("the rest time: " + str(t3 - t2))
            return str(netcdf_file_folder)
                #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))

        # then add the network link using this args.server
        elif len(endpoint_list) == 0:
            #empty_folder = kml.newfolder(name="no points in view")
            #return str(empty_folder)
            pass
        # for surveys with 1 or 2 points. Can't make a polygon.
        else:
            print("not enough points")
            #empty_folder = kml.newfolder(name="no points in view")
            #return str(empty_folder)
            # single
            #print("single survey")


    else:  # polygons only

        if len(endpoint_list) > 0:
             netcdf_file_folder = kml.newfolder()
             for netcdf in endpoint_list:
                print("NETCDF: " + str(netcdf))

                converter_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
                print(converter_obj.npu.point_count)
                if converter_obj.npu.point_count > 2:
                    polygon_folder, polygon = converter_obj.build_polygon(netcdf_file_folder)
                    dataset_polygon_region = converter_obj.build_region(converter_obj.MAX_LOD_PIXELS, -1,
                                                               converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)
                    polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder

                else:  # for surveys with 1 or 2 points. Can't make a polygon. Still save the points?
                    print("not enough points")
                    #empty_folder = kml.newfolder(name="no points in view")
                    #return str(empty_folder)
             return str(netcdf_file_folder)

        else:
            empty_folder = kml.newfolder(name="no points in view")
            return str(empty_folder)


app.run(debug=True)

