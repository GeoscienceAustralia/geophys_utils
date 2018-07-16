from flask import Flask
from flask_restful import Api, Resource, reqparse
from flask import request
import simplekml
from geophys_utils import NetCDFPointUtils, get_spatial_ref_from_wkt
from geophys_utils import points2convex_hull
from geophys_utils import _polygon_utils
import netCDF4


app = Flask(__name__)
api = Api(app)

@app.route('/<bounding_box>',methods=['GET'])

def generate_subset_kml_points(bounding_box):

    query_string = request.args
    bbox = query_string['BBOX']
    bbox_list = bbox.split(',')
    west = float(bbox_list[0])
    south = float(bbox_list[1])
    east = float(bbox_list[2])
    north = float(bbox_list[3])
    new_west = ((east - west) / 2) + west
    print(new_west)
    # Create NetCDFPointUtils object for specified netCDF dataset

    netcdf_path = "http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc"
    print(netcdf_path)
    # netcdf_path = 'E:\\Temp\\gravity_point_test\\195256.nc'
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

    spatial_mask = npu.get_spatial_mask([new_west, south, east, north])
    #lookup_var = npu.get_lookup_mask(['obsno', 'grav'])
    #print(lookup_var)
    if True in spatial_mask:
        # when ordered through the all_point_data_generator it appears to come out as [obsno, lat, long, then everything else as in field_list]
        field_list = ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype', 'reliab']#, 'freeair', '', '']
        print(spatial_mask)
        point_data_generator = npu.all_point_data_generator(field_list, spatial_mask)
        print(len(spatial_mask))

        variable_attributes = next(point_data_generator)
        # print(variable_attributes.keys())
        # print(variable_attributes)

        # Use long names instead of variable names where they exist
        field_names = [variable_attributes[variable_name].get('long_name') or variable_name for variable_name in
                       variable_attributes.keys()]
        print('field_names: {}'.format(field_names))

        skip_points = 1
        points_folder = kml.newfolder()
        points_read = 0
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
            print("POINT DATA")
            print(point_data)
            new_point = points_folder.newpoint(name=str(point_data[0]), coords=[(point_data[2], point_data[1])])

            description_string = '<![CDATA[' \
                                 '<h3>Point Data:</h3>' \
                                 '<p><b>{0}: </b>{1}</p>' \
                                 '<p><b>{2}: </b>{3}</p> ' \
                                 '<p><b>{4}: </b>{5}</p>' \
                                 '<p><b>{6}: </b>{7}</p> ' \
                                 '<p><b>{8}: </b>{9}</p>' \
                                 ']]>'.format(field_names[3], point_data[3],
                                              field_names[4], point_data[4],
                                              field_names[5], point_data[5],
                                              field_names[6], point_data[6],
                                              field_names[7], point_data[7]
                                              )
            print(description_string)
            new_point.description = description_string

            # set the point icon. Different urls can be found in point style options in google earth
            new_point.style.iconstyle.icon.href = POINT_ICON_STYLE_LINK
            new_point.style.iconstyle.scale = 0.7
            new_point.labelstyle.scale = 0  # removes the label


        print(points_folder)

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
