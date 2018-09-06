# ===============================================================================
#    Copyright 2017 Geoscience Australia
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
# ===============================================================================
'''
Created on 03/08/2018

@author: Andrew Turner
'''
import simplekml
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging
from datetime import date
from geophys_utils.dataset_metadata_cache import SQLiteDatasetMetadataCache
import netCDF4
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils
import yaml
import os
import numpy as np

# TODO: GET RID OF THIS HORRIBLE HACK
# Set the following to None or empty string to use OPeNDAP endpoints
# LOCAL_FILE_LOCATION = None
LOCAL_FILE_LOCATION = 'C:\\Users\\u62231\\Desktop\\aem_data\\AusAEM_Year1_Tranche1_Final_EM.nc'
# LOCAL_FILE_LOCATION = 'C:\\Users\\u62231\\Desktop\\grav_data_10_july'

COLORMAP_NAME = 'rainbow'
COLOUR_STRETCH_RANGE = (-500, 500)  # min/max tuple for colour stretch range

# Setup logging handlers if required
logger = logging.getLogger(__name__)  # Get logger
logger.setLevel(logging.DEBUG)  # Initial logging level for this module

LINE_RESOLUTION_STEPS = 70


def convert_value_from_old_to_new_range(value_to_convert, old_range_min, old_range_max, new_range_min, new_range_max):
    # converts a value from a np array to a value within a desired range. Essentially it converts a number in one
    # range to a number in another range, while maintaining the ratio.
    old_min = old_range_min
    old_range = old_range_max - old_range_min
    new_range = new_range_max - new_range_min

    new_value = (((value_to_convert - old_min) * new_range) / old_range) + 0

    return new_value


class NetCDF2kmlConverter(object):
    def __init__(self, netcdf_path, dataset_settings, metadata_tuple=None):
        logger.debug(metadata_tuple)
        self.npu = None
        self.survey_id = metadata_tuple[0]
        self.survey_title = metadata_tuple[1]
        self.netcdf_path = metadata_tuple[2]  # test
        self.polygon = metadata_tuple[3]

        self.kml = simplekml.Kml()

        # Store dataset spatial extents as python variables
        self.west_extent = metadata_tuple[4]
        self.east_extent = metadata_tuple[5]
        self.south_extent = metadata_tuple[6]
        self.north_extent = metadata_tuple[7]

        self.point_count = metadata_tuple[8]

        self.start_date = metadata_tuple[9]
        self.end_date = metadata_tuple[10]
        self.colormap = plt.cm.get_cmap(COLORMAP_NAME, 256)

        # TODO: GET RID OF THIS HORRIBLE HACK
        # if LOCAL_FILE_LOCATION:
        #     self.netcdf_path = os.path.join(LOCAL_FILE_LOCATION,
        #                            os.path.basename(self.netcdf_path)
        #                            )

        self.netcdf_path = netcdf_path

        self.netcdf_dataset = netCDF4.Dataset(self.netcdf_path)
        self.npu = NetCDFPointUtils(self.netcdf_dataset)
        self.thredds_metadata_link = dataset_settings['thredds_link']

        self.polygon_style = simplekml.Style()
        self.polygon_style.polystyle.color = dataset_settings['polygon_color']
        self.polygon_style.polystyle.outline = dataset_settings['polygon_outline']

        self.point_flag_exists = dataset_settings['point_flag_exists']

        self.point_style = simplekml.Style()
        self.point_style.iconstyle.scale = dataset_settings['point_icon_scale']
        self.point_style.labelstyle.scale = dataset_settings['point_labelstyle_scale']  # 0 removes the label
        self.point_style.iconstyle.icon.href = dataset_settings['point_icon_href']
        self.unflagged_point_icon_color_scheme = None  # insert function string
        self.flagged_point_icon_color = 'ff000000'

        self.field_list = dataset_settings['point_field_list']

    def build_region(self, min_lod_pixels=100, max_lod_pixels=-1, min_fade_extent=200, max_fade_extent=800):
        """
        Builds a KML Region using simplekml.
        lod parameters are the measurements in screen pixels that represents the maximum limit of the visibility range
        for a given Region.
        :param min_lod_pixels:
        :param max_lod_pixels: -1 the default, indicates "active to infinite size."
        :param min_fade_extent:
        :param max_fade_extent:
        :return: region object in simplekml
        """

        region = simplekml.Region(latlonaltbox="<north>" + str(self.north_extent) + "</north>" +
                                               "<south>" + str(self.south_extent) + "</south>" +
                                               "<east>" + str(self.east_extent) + "</east>" +
                                               "<west>" + str(self.west_extent) + "</west>",
                                  lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
                                      "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
                                      "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
                                      "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")
        return region

    def build_polygon(self, parent_folder, visibility=True):
        """
        Builds a kml polygon into the parent folder. Polygon is build from netcdf flobal attribute geospatial_bounds.
        :param parent_folder:
        :param visibility:
        :return: the input parent folder now containing the polygon and the polygon itself if that is desired instead.
        """
        logger.debug("Building polygon...")
        # get the polygon points from the netcdf global attributes
        try:
            polygon_bounds = [[float(ordinate)
                               for ordinate in coord_pair.strip().split(' ')
                               ]
                              for coord_pair in
                              re.search('POLYGON\(\((.*)\)\)',
                                        self.polygon
                                        ).group(1).split(',')
                              ]
            # build the polygon based on the bounds. Also set the polygon name. It is inserted into the parent_folder.
            pol = parent_folder.newpolygon(name=str(self.survey_title) + " " + str(self.survey_id),
                                           outerboundaryis=polygon_bounds, visibility=visibility)

            # build the polygon description
            description_string = '<![CDATA['
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name',
                                                                                      str(self.survey_title))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', str(self.survey_id))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Start Date',
                                                                                      str(self.start_date))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey End Date',
                                                                                      str(self.end_date))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('NCI Data Link', str(
                self.thredds_metadata_link) + str(self.survey_id) + '.nc')
            description_string = description_string + ']]>'
            pol.description = description_string

            # Add timestamp
            # if self.start_date is not None and self.end_date is not None:
            try:
                assert self.start_date > date(1900, 1, 1), 'Start date {} is less than 1900-01-01'.format(
                    self.start_date)
                assert self.end_date < date(2020, 1, 1), 'End date {} is greater than 2020-01-01'.format(self.end_date)
                assert self.end_date > date(1900, 1, 1), 'End date {} is less than 1900-01-01'.format(self.end_date)
                assert self.start_date < date(2020, 1, 1), 'Start date {} is greater than 2020-01-01'.format(
                    self.start_date)
                pol.timespan.begin = str(self.start_date)
                pol.timespan.end = str(self.end_date)

            except:  # if survey does not contain start/end date information, use survey id as start/end date year.
                survey_year = int(re.match('^[0-9]{4}', str(self.survey_id)).group())
                assert survey_year > 1900 and survey_year < 2020, 'survey_year <= 1900 or survey_year >= 2020'

                pol.timespan.begin = str(survey_year) + "-06-01"
                pol.timespan.end = str(survey_year) + "-07-01"
            pol.style = self.polygon_style

        except Exception as e:
            logger.warning("Unable to display polygon: {}".format(e))

        return parent_folder

    def build_lines(self, netcdf_file_folder, bounding_box):

        logger.debug("Building lines...")
        bounding_box_floats = [float(coord) for coord in bounding_box]
        line_utils = NetCDFLineUtils(self.netcdf_dataset)
        line_list = [line for line in self.netcdf_dataset.variables['line']]
        logger.debug("List of lines in view: {}".format(line_list))

        line_dict = line_utils.get_lines(line_list, variables=['gps_elevation', "line_index"], bounds=bounding_box_floats)
        logger.debug("Line Dictionary: {}".format(next(line_dict)))

        for line in line_dict:
            line_number = line[0]
            line_folder = netcdf_file_folder.newfolder(name="Line Number: {}".format(line_number))
            coords_array = line[1]['coordinates']
            #number_of_points_in_line = int(coords_array.size / 2)
            number_of_points_in_line = len(line[1]['gps_elevation'])

            if coords_array.size is not 0:

                small_array = coords_array[0:coords_array.size:LINE_RESOLUTION_STEPS]  # step size is hardcoded at 100

                # build the line segments
                i = 1
                elevation_index = LINE_RESOLUTION_STEPS
                while i < small_array.size / 2:
                    last_point_index = int(small_array.size)
                    linestring = line_folder.newlinestring(name=str("Line Segment: " + str(i - 1)))

                    linestring.coords = [
                                         (small_array[i - 1][0], small_array[i - 1][1], line[1]['gps_elevation'][elevation_index - LINE_RESOLUTION_STEPS]),
                                         (small_array[i][0], small_array[i][1], line[1]['gps_elevation'][elevation_index])
                                        ]
                    linestring.altitudemode = simplekml.AltitudeMode.relativetoground
                    #linestring.altitudemode = simplekml.AltitudeMode.clamptoground
                    linestring.extrude = 0
                    linestring.tessellate = 1
                    linestring.style.linestyle.width = 5
                    linestring.style.linestyle.color = simplekml.Color.yellowgreen

                    # build the polygon description
                    description_string = '<![CDATA['
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Line Number', str(line_number))
                    description_string = description_string + ']]>'
                    linestring.description = description_string

                    i = i + 1
                    elevation_index = elevation_index + LINE_RESOLUTION_STEPS

                #linestring = lines_folder.newlinestring(name="Line Number: {}".format(line_number))
                linestring.coords = [(small_array[i - 1][0], small_array[i - 1][1], line[1]['gps_elevation'][elevation_index - LINE_RESOLUTION_STEPS]),
                                     (coords_array[-1][0], coords_array[-1][1], line[1]['gps_elevation'][number_of_points_in_line- 1])
                                     ]
                # # # touring
                # tour = kml.newgxtour(name="Play me!")
                # playlist = tour.newgxplaylist()
                #
                # soundcue = playlist.newgxsoundcue()
                # soundcue.href = "http://code.google.com/p/simplekml/source/browse/samples/drum_roll_1.wav"
                # soundcue.gxdelayedstart = 2
                #
                # animatedupdate = playlist.newgxanimatedupdate(gxduration=6.5)
                # animatedupdate.update.change = '<IconStyle targetId="{0}"><scale>10.0</scale></IconStyle>'.format(
                #     pnt.style.iconstyle.id)
                #
                # flyto = playlist.newgxflyto(gxduration=4.1)
                # flyto.camera.longitude = 170.157
                # flyto.camera.latitude = -43.671
                # flyto.camera.altitude = 9700
                # flyto.camera.heading = -6.333
                # flyto.camera.tilt = 33.5
                # flyto.camera.roll = 0

                # wait = playlist.newgxwait(gxduration=2.4)

            else:
                print("line doesn't have any points in view")

        return netcdf_file_folder

    def build_points(self, points_folder, bounding_box):
        """
        Builds all points for a survey. Including building the containing folder, setting time stamps, setting the
         style, and setting the description html to pop up when the point is selected.
        :param points_folder: The folder for the new survey folder containing the points to be inserted into.
        :return: the kml folder of the survey containing the new points.
        """
        logger.debug("Building points for netcdf file: " + str(self.netcdf_path))
        logger.debug('bounding_box:' + str(bounding_box))
        bounding_box_floats = [float(coord) for coord in bounding_box]
        print(bounding_box_floats)
        # bounding_box_floats = [110.4899599829594, -56.11642075733719, 166.658146968822, -6.11642075733719]
        spatial_mask = self.npu.get_spatial_mask(bounding_box_floats)
        print(spatial_mask)
        if True in spatial_mask:
            print("TRUE")
            new_survey_folder = points_folder.newfolder(name=str(self.survey_title) + " " + str(self.survey_id))

            # Set timestamp
            # start_date = re.match('^[0-9]{4}', str(self.survey_id)).group()
            # new_survey_folder.timespan.begin = str(start_date) + "-01-01"
            # new_survey_folder.timespan.end = str(start_date) + "-01-01"

            point_data_generator = self.npu.all_point_data_generator(self.field_list, spatial_mask)
            logger.debug(point_data_generator)
            variable_attributes = next(point_data_generator)
            logger.debug(variable_attributes)
            print("variable_attributes: " + str(variable_attributes))
            aem_variable_atts = {'aem_test': {'units': 'um/s^2'}}

            skip_points = 1  # set to limit the points displayed if required.
            points_read = 0

            for point_data in point_data_generator:
                print("POINT DATA: " + str(point_data))
                points_read += 1

                # ignore points between skip_points
                if points_read % skip_points != 0:
                    continue

                # add new points with netcdf file Obsno as title and long and lat as coordinatess
                new_point = new_survey_folder.newpoint(name="Observation no. " + str(point_data[0]),
                                                       coords=[(point_data[2], point_data[1])])

                description_string = self.build_html_description_string(variable_attributes, point_data)
                logger.debug(description_string)
                new_point.description = description_string  # set description to point

                # styling
                self.unflagged_point_icon_color_scheme = self.value2colourhex(point_data[5])
                new_point.style = self.point_style  # sets the icon scale, labelstyle scale, icon link

                # Set the style for points if gridflag == 2
                if self.point_flag_exists:  # if there is a point_flag separate the points and colour differently
                    if point_data[8] == "Station not used in the production of GA grids.":
                        new_point.style.iconstyle.color = self.flagged_point_icon_color

            return new_survey_folder
        else:
            return None

    def build_html_description_string(self, variable_attributes, point_data):
        """
            Helper function to build the description string automatically based on how many fields are in the field list.
            It will take into account if a unit of measure needs to be specified or not.
            :param field_list: The fields to be included in the description.
            :param variable_attributes: a list of lists? Matching point_data.
            :param point_data: A list of the actual data of each point
            :return: return the description string html ready to be attached to a kml point.
        """

        print(variable_attributes)

        description_string = '<![CDATA['
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', self.survey_title)
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', self.survey_id)

        i = 1  # skip obsno, lat, long in field_list
        print("LENGTH: " + str(len(self.field_list)))
        while i < len(self.field_list) - 1:
            print(self.field_list)
            print(point_data[i])
            if variable_attributes[self.field_list[i]].get('units'):
                description_string = description_string + '<p><b>{0}: </b>{1} {2}</p>'.format(
                    variable_attributes[self.field_list[i]].get('long_name'),
                    point_data[i],
                    variable_attributes[self.field_list[i]].get('units'))
            else:
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format(
                    variable_attributes[self.field_list[i]].get('long_name'),
                    point_data[i])
            i += 1
            print(i)
        description_string = description_string + ']]>'
        return description_string

    def get_basic_density_measurement(self):
        """
        Basic calculation to get an idea of the dataset point density. This could be used to dynamically set the
        min_lod_pixels, max_lod_pixels, and point scale to get the best visual outcome. A higher density may need a smaller
        point scale or even an additional region.
        :return: Currently just prints to output
        """
        logger.debug("Point Count: {}".format(self.npu.point_count))
        logger.debug("Area: {}".format((self.east_extent - self.west_extent) * (self.north_extent - self.south_extent)))
        logger.debug("Density: {}".format(((self.east_extent - self.west_extent) * (
        self.north_extent - self.south_extent)) / self.npu.point_count * 1000))

    def build_static_kml(self):

        point_style = simplekml.Style()
        point_style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
        point_style.iconstyle.scale = 0.7
        point_style.labelstyle.scale = 0  # removes the label

        bbox = ['115.6900403947253', '-80.3741299166767', '165.9763920890655', '-7.21307533348337']
        # polygon style
        polygon_style = simplekml.Style()
        polygon_style.polystyle.color = '990000ff'  # Transparent red
        polygon_style.polystyle.outline = 1

        kml = simplekml.Kml()
        netcdf_file_folder = kml.newfolder(name=self.survey_title + " " + self.survey_id)

        polygon_folder = self.kml.newfolder(name="polygon")
        polygon_folder, polygon = self.build_polygon(polygon_folder)
        dataset_polygon_region = self.build_region(-1, 600, 200, 800)
        dataset_points_region = self.build_region(0, -1, 200, 800)

        points_folder = netcdf_file_folder.newfolder(name="points")
        points_folder = self.build_points(points_folder, bbox)

        # structure them correctly
        polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
        points_folder.region = dataset_points_region  # insert built point region into point folder
        return self.kml

    def value2colourhex(self, data_value):
        '''
        Function to convert data value to hex string for color in self.colormap
        '''
        cmap_value = int(convert_value_from_old_to_new_range(data_value, *COLOUR_STRETCH_RANGE, 0, 255))
        # print(cmap_value)
        # print(str(self.colormap(cmap_value)[:3]))
        # print(str(mpl.colors.rgb2hex(self.colormap(cmap_value)[:3])))
        # print(mpl.colors.rgb2hex(self.colormap(cmap_value)[:3]))

        return mpl.colors.rgb2hex(self.colormap(cmap_value)[:3]).replace('#', 'ff')


def build_dynamic_network_link(containing_folder, link="http://127.0.0.1:5000/query"):
    """
    Build a network link, set the parameters, and inserts into the specified containing folder.
    """
    net_link = containing_folder.newnetworklink(name="Network Link")
    net_link.link.href = link
    net_link.link.viewrefreshmode = simplekml.ViewRefreshMode.onstop
    net_link.link.viewrefreshtime = 1
    net_link.link.refreshinterval = 2

    return net_link


def main():
    # --------------------------------
    #   build dynamic_grav kml
    # --------------------------------
    # kml = simplekml.Kml()
    # new_folder = kml.newfolder()
    # build_dynamic_network_link(new_folder)
    # kml.save("dynamic_grav_surveys.kml")

    # --------------------------------
    # build static kml
    # --------------------------------
    bounding_box = [[100.00, -50.00], [159.00, -5.00]]  # include all of aus
    sdmc = SQLiteDatasetMetadataCache(debug=False)
    endpoint_list = sdmc.search_dataset_distributions(
        keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
        protocol='opendap',
        ll_ur_coords=bounding_box
    )
    for survey in endpoint_list:
        kml = simplekml.Kml()
        netcdf_file_folder = kml.newfolder(name="Ground Gravity Survey Observations")
        netcdf2kml_obj = NetCDF2kmlConverter(survey)
        netcdf2kml_obj.netcdf_dataset = netCDF4.Dataset(netcdf2kml_obj.netcdf_path)
        netcdf2kml_obj.npu = NetCDFPointUtils(netcdf2kml_obj.netcdf_dataset)

        survey_region = netcdf2kml_obj.build_region()
        polygon_region = netcdf2kml_obj.build_region(0, 100)

        polygon_folder = kml.newfolder(name="polygon")
        polygon_folder = netcdf2kml_obj.build_polygon(polygon_folder)
        polygon_folder.region = polygon_region  # insert built point region into point folder

        netcdf_file_folder.region = survey_region  # insert built point region into point folder

        if netcdf2kml_obj.npu.point_count > 2:

            survey_points_folder = netcdf2kml_obj.build_points(netcdf_file_folder,
                                                               ['110.00', '-45.00', '155.00', '-10.00'])
            survey_points_folder.region = survey_region  # insert built point region into point folder

            cleaned_survey_title = re.sub("/", " or ", netcdf2kml_obj.survey_title)
            cleaned_survey_title = re.sub('"', "", cleaned_survey_title)

            kml.save('C:\\Users\\u62231\\Desktop\\grav_kmls\\' + str(cleaned_survey_title) + ' ' +
                     str(netcdf2kml_obj.survey_id) + '.kml')
        else:
            print("fail")
            netcdf2kml_obj.netcdf_dataset.close()



            # netcdf_path_list = []
            # i = 2
            # while i < len(sys.argv):
            #     print(sys.argv[i])
            #     netcdf_path_list.append(sys.argv[i])
            #     i += 1
            # print(netcdf_path_list)

            # parser = argparse.ArgumentParser()
            #
            # parser.add_argument("-s", "--server", help="The server to receive the get request from google earth and dynamically "
            #                                            "build kml points within the bbox. If this parameter is empty, a static "
            #                                            "kml will be generated", type=str, required=False)
            # parser.add_argument("-n", "--netcdf_path_list", help="Add one or more paths to netcdf files to be converted into a"
            #                                                      "single kml file.", type=str, nargs='+')
            #
            # args = parser.parse_args()
            #
            # print(args.server)
            # print(args.netcdf_path_list)

            # if args.server:
            #     # dynamic build
            #     print("dynamic build")
            #     if len(args.netcdf_path_list) > 1:
            #         print("multiple surveys")
            #         # multiples
            #         list_of_converter_objects= []
            #         for netcdf in args.netcdf_path_list:
            #             #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))
            #             pass
            #
            #         # then add the network link using this args.server

            #     else:
            #         # single
            #         print("single survey")
            #         converter_object = NetCDF2kmlConverter(args.netcdf_path_list[0])
            #         converter_object.build_dynamic_kml()
            #         converter_object.kml.save(converter_object.survey_title + " dynamic points.kml")
            #         print("Building kml for survey: " + converter_object.survey_title + " dynamic points.kml")
            # else:
            #     # static build
            #     print("static build")
            #     if len(args.netcdf_path_list) > 1:
            #         print("multiple surveys")
            #         pass
            #     else:
            #         print("single survey")
            #         converter_object = NetCDF2kmlConverter(args.netcdf_path_list[0])
            #         converter_object.build_static_kml()
            #
            #         converter_object.kml.save(converter_object.survey_title + " static points.kml")
            #         print("Building kml for survey: " + converter_object.survey_title + " static points.kml")
            # single

            # modified_server_url_for_dynamic_kml_generation = args.server + 'query'
            # NetCDF2kmlConverter(args.netcdf_path_list)


if __name__ == '__main__':
    main()
