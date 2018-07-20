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

app = Flask(__name__)
api = Api(app)


def do_the_things(netcdf_link, points_folder):
    t0 = time.time()  # retrieve coordinates from query
    query_string = request.args
    print('query done!')
    print(time.time())
    bbox = query_string['BBOX']
    bbox_list = bbox.split(',')
    west = float(bbox_list[0])
    south = float(bbox_list[1])
    east = float(bbox_list[2])
    north = float(bbox_list[3])
    print(time.time())

    t1 = time.time()  # Create NetCDFPointUtils object for specified netCDF dataset
    # netcdf_path = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc"
    # 195256
    # netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_data_10_july\\201780.nc"
    # 195105
    # 201780

    POINT_ICON_STYLE_LINK = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
    netcdf_dataset = netCDF4.Dataset(netcdf_link)
    npu = NetCDFPointUtils(netcdf_dataset)

    survey_title = str(npu.netcdf_dataset.getncattr('title'))
    kml = simplekml.Kml()

    # REGIONS



    # POLYGONS



    # POINTS

    t2 = time.time()  # create the spatial mask
    spatial_mask = npu.get_spatial_mask([west, south, east, north])

    t3 = time.time()  # get the points and variable info from point generator
    if True in spatial_mask:
        # when ordered through the all_point_data_generator it appears to come out as [obsno, lat, long, then everything else as in field_list]
        field_list = ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype',
                      'reliab']  # , 'freeair', '', '']

        point_data_generator = npu.all_point_data_generator(field_list, spatial_mask)
        variable_attributes = next(point_data_generator)

        # Use long names instead of variable names where they exist
        # field_names = [variable_attributes[variable_name].get('long_name') or variable_name for variable_name in
        #                variable_attributes.keys()]

        skip_points = 3
        points_read = 0
        t4 = time.time()  # loop through the points and create them.
        for point_data in point_data_generator:
            points_read += 1

            # ignore points between skip_points
            if points_read % skip_points != 0:
                continue

            # <b>Observation Number: </b><br>
            # key long name: value <br>
            #

            # print(point_data[0])
            # add new points with netcdf file Obsno as title and long and lat as coordinatess

            new_point = points_folder.newpoint(name="Point no. " + str(point_data[0]),
                                               coords=[(point_data[2], point_data[1])])

            description_string = '<![CDATA[' \
                                 '<p><b>{0}: </b>{1} {2}</p>' \
                                 '<p><b>{3}: </b>{4} {5}</p> ' \
                                 '<p><b>{6}: </b>{7} {8}</p>' \
                                 '<p><b>{9}: </b>{10}</p> ' \
                                 '<p><b>{11}: </b>{12}</p>' \
                                 ']]>'.format(
                variable_attributes['grav'].get('long_name'), point_data[3], variable_attributes['grav'].get('units'),
                variable_attributes['freeair'].get('long_name'), point_data[4],
                variable_attributes['freeair'].get('units'),  # free air
                variable_attributes['bouguer'].get('long_name'), point_data[5],
                variable_attributes['bouguer'].get('units'),  # bouguer
                variable_attributes['stattype'].get('long_name'), point_data[6],  # station type
                variable_attributes['reliab'].get('long_name'), point_data[7]  # reliability
            )

            new_point.description = description_string

            # set the point icon. Different urls can be found in point style options in google earth
            new_point.style.iconstyle.icon.href = POINT_ICON_STYLE_LINK
            new_point.style.iconstyle.scale = 0.7
            new_point.labelstyle.scale = 0  # removes the label

        t5 = time.time()

        time_get_query_points = t1 - t0
        time_create_netcdf_object = t2 - t1
        time_create_spatial_mask = t3 - t2
        time_point_gen = t4 - t3
        time_create_points = t5 - t4

        print("time_get_query_points: " + str(time_get_query_points))
        print("time_create_netcdf_object: " + str(time_create_netcdf_object))
        print("time_create_spatial_mask: " + str(time_create_spatial_mask))
        print("time_point_gen: " + str(time_point_gen))
        print("time_create_points: " + str(time_create_points))
        print(points_folder)
        return points_folder
    else:
        print("no points in view")

        return "<Folder><name>No points in view</name></Folder>"


@app.route('/<bounding_box>',methods=['GET'])
def do_everything(bounding_box):

    query_string = request.args
    parser = argparse.ArgumentParser()

    # parser.add_argument("-s", "--server",
    #                     help="The server to receive the get request from google earth and dynamically "
    #                          "build kml points within the bbox. If this parameter is empty, a static "
    #                          "kml will be generated", type=str, required=False)
    parser.add_argument("-n", "--netcdf_path_list", help="Add one or more paths to netcdf files to be converted into a"
                                                         "single kml file.", type=str, nargs='+')

    args = parser.parse_args()

    kml = simplekml.Kml()
    points_folder = kml.newfolder()

    if len(args.netcdf_path_list) > 1:
        print("multiple surveys")
        # multiples
        list_of_point_folders= []



        for netcdf in args.netcdf_path_list:
            netcdf_file_folder = kml.newfolder()

            converter_obj = netcdf2kml.NetCDF2kmlConverter(netcdf)
            # region = converter_obj.build_region()
            polygon_folder, polygon = converter_obj.build_polygon(netcdf_file_folder)

            # dataset_points_folder = plot_points_in_kml()
            # dataset_points_folder = plot_dynamic_points_in_kml()
            dataset_polygon_region = converter_obj.build_region(converter_obj.MAX_LOD_PIXELS, converter_obj.MIN_LOD_PIXELS + 400,
                                                       converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)
            dataset_points_region = converter_obj.build_region(converter_obj.MIN_LOD_PIXELS - 200, converter_obj.MAX_LOD_PIXELS,
                                                               converter_obj.MIN_FADE_EXTENT, converter_obj.MAX_FADE_EXTENT)
            # polygon_folder, polygon = converter_obj.build_polygon()
            #converter_obj.build_static_kml()
            polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
            points_folder.region = dataset_points_region  # insert built point region into point folder
            points_folder = converter_obj.do_the_things(netcdf, points_folder, query_string)
            #points_folder = do_the_things(netcdf, points_folder)
        kml.save("test.kml")
        return str(points_folder)
            #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))

    # then add the network link using this args.server

    else:
        # single
        print("single survey")

        folder = do_the_things(args.netcdf_path_list[0], points_folder)
        # converter_object = netcdf2kml.NetCDF2kmlConverter(args.netcdf_path_list[0])
        # converter_object.build_dynamic_kml()
        # converter_object.kml.save(converter_object.survey_title + " dynamic points.kml")
        # print("Building kml for survey: " + converter_object.survey_title + " dynamic points.kml")
        return str(folder)

app.run(debug=True)

