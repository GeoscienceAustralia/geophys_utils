from flask import Flask
from flask_restful import Api, Resource, reqparse
from flask import request
import simplekml
from geophys_utils import NetCDFPointUtils, get_spatial_ref_from_wkt
from geophys_utils import points2convex_hull
from geophys_utils import _polygon_utils
import netCDF4
import time


app = Flask(__name__)
api = Api(app)

@app.route('/<bounding_box>',methods=['GET'])

def generate_subset_kml_points(bounding_box):

    t0 = time.clock()  # retrieve coordinates from query
    query_string = request.args
    bbox = query_string['BBOX']
    bbox_list = bbox.split(',')
    west = float(bbox_list[0])
    south = float(bbox_list[1])
    east = float(bbox_list[2])
    north = float(bbox_list[3])
    new_west = ((east - west) / 2) + west

    t1 = time.clock()     # Create NetCDFPointUtils object for specified netCDF dataset
    netcdf_path = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc"
    # 195256
    # netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_netcdf_4\\195107.nc"
    # 195105
    # 201780

    POINT_ICON_STYLE_LINK = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
    netcdf_dataset = netCDF4.Dataset(netcdf_path)
    npu = NetCDFPointUtils(netcdf_dataset)

    survey_title = str(npu.netcdf_dataset.getncattr('title'))
    kml = simplekml.Kml()
# -----------------------bounds------------------
#   self.bounds = [xmin, ymin, xmax, ymax]


    t2 = time.clock()  # create the spatial mask
    spatial_mask = npu.get_spatial_mask([new_west, south, east, north])

    t3 = time.clock() # get the points and variable info from point generator
    if True in spatial_mask:
        # when ordered through the all_point_data_generator it appears to come out as [obsno, lat, long, then everything else as in field_list]
        field_list = ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype', 'reliab']#, 'freeair', '', '']

        point_data_generator = npu.all_point_data_generator(field_list, spatial_mask)
        variable_attributes = next(point_data_generator)


        # Use long names instead of variable names where they exist
        # field_names = [variable_attributes[variable_name].get('long_name') or variable_name for variable_name in
        #                variable_attributes.keys()]



        points_folder = kml.newfolder()

        skip_points = 1
        points_read = 0
        t4 = time.clock() # loop through the points and create them.
        for point_data in point_data_generator:
            points_read += 1

            # ignore points between skip_points
            if points_read % skip_points != 0:
                continue

            #<b>Observation Number: </b><br>
            # key long name: value <br>
            #

            #print(point_data[0])
            # add new points with netcdf file Obsno as title and long and lat as coordinatess


            new_point = points_folder.newpoint(name="Point no. " + str(point_data[0]), coords=[(point_data[2], point_data[1])])

            description_string = '<![CDATA[' \
                                 '<p><b>{0}: </b>{1} {2}</p>' \
                                 '<p><b>{3}: </b>{4} {5}</p> ' \
                                 '<p><b>{6}: </b>{7} {8}</p>' \
                                 '<p><b>{9}: </b>{10}</p> ' \
                                 '<p><b>{11}: </b>{12}</p>' \
                                 ']]>'.format(
                variable_attributes['grav'].get('long_name'), point_data[3], variable_attributes['grav'].get('units'),
                variable_attributes['freeair'].get('long_name'), point_data[4], variable_attributes['freeair'].get('units'),  # free air
                variable_attributes['bouguer'].get('long_name'), point_data[5], variable_attributes['bouguer'].get('units'),  # bouguer
                variable_attributes['stattype'].get('long_name'), point_data[6],  # station type
                variable_attributes['reliab'].get('long_name'), point_data[7]  # reliability
                                              )

            new_point.description = description_string

            # set the point icon. Different urls can be found in point style options in google earth
            new_point.style.iconstyle.icon.href = POINT_ICON_STYLE_LINK
            new_point.style.iconstyle.scale = 0.7
            new_point.labelstyle.scale = 0  # removes the label


        print(points_folder)

        t5 = time.clock()

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

        return str(points_folder)
    else:
        print("no points in view")

        return "<Folder><name>No points in view</name></Folder>"




    # def getBbox(bounding_box):
    #     #if request.method == "GET":
    #         #if re.search(bounding_box, "query"):
    #
    #             print('si')
    #             print(request.query_string)
    #             query_string = request.args
    #             #request.
    #             # print(bounding_box)
    #             # #t = re.search('=*.', bounding_box)
    #             # t = re.split("=", bounding_box)
    #             # print(t[1])
    #             #print(query_string['BBOX'])
    #             bbox = query_string['BBOX']
    #             bbox_list = bbox.split(',')
    #             west = float(bbox_list[0])
    #             south = float(bbox_list[1])
    #             east = float(bbox_list[2])
    #             north = float(bbox_list[3])
    #             #print("north: " + str(north))
    #
    #             center_lng = ((east - west) / 2) + west
    #             center_lat = ((north - south) / 2) + south
    #
    #             kml =     '<?xml version="1.0" encoding="UTF-8"?>' \
    #                       '<kml xmlns="http://www.opengis.net/kml/2.2">' \
    #                       '<Placemark>' \
    #                       '<name>View-centered placemark</name>' \
    #                       '<Point>' \
    #                       '<coordinates>{0}.6f,{1}.6f</coordinates>' \
    #                       '</Point>' \
    #                       '</Placemark>' \
    #                       '</kml>'.format(center_lng, center_lat)
    #             print('Content-Type: application/vnd.google-earth.kml+xml\n')
    #             print(kml)
    #             #return "<Folder><name>" + kml + "</name></Folder>"
    #             #return "application/vnd.google-earth.kml+xml"
    #             #return 'Content-Type: application/vnd.google-earth.kml+xml\n', kml
    #             return kml




app.run(debug=True)
