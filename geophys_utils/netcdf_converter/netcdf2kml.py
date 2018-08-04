import simplekml
import re
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import logging

COLORMAP_NAME='rainbow'
COLOUR_STRETCH_RANGE=(-500, 500) # min/max tuple for colour stretch range

# Setup logging handlers if required
logger = logging.getLogger(__name__)  # Get logger
logger.setLevel(logging.INFO)  # Initial logging level for this module


def convert_value_from_old_to_new_range(value_to_convert, old_range_min, old_range_max, new_range_min, new_range_max):
# converts a value from a np array to a value within a desired range. Essentially it converts a number in one
# range to a number in another range, while maintaining the ratio.
    old_min = old_range_min
    old_range = old_range_max - old_range_min
    new_range = new_range_max - new_range_min

    new_value = (((value_to_convert - old_min) * new_range) / old_range) + 0

    return new_value

class NetCDF2kmlConverter(object):

    def __init__(self, metadata_tuple=None):


        # Create NetCDFPointUtils object for specified netCDF dataset
        #netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc'
        #netcdf_path = 'E:\\Temp\\gravity_point_test\\195256.nc'
        #195256
        #netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_netcdf_4\\201780.nc"
        #195105
        #201780

        self.npu = None
        self.survey_id = metadata_tuple[0]
        self.survey_title = metadata_tuple[1]
        self.netcdf_path = metadata_tuple[2]
        self.polygon = metadata_tuple[3]

        self.kml = simplekml.Kml()

        # Store dataset spatial extents as python variables
        self.west_extent = metadata_tuple[4]
        self.east_extent = metadata_tuple[5]
        self.south_extent = metadata_tuple[6]
        self.north_extent = metadata_tuple[7]

        self.point_icon_style_link = "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"
        
        self.colormap = plt.cm.get_cmap(COLORMAP_NAME, 256)

    # -----------------------bounds------------
    #        self.bounds = [xmin, ymin, xmax, ymax]

    def build_region(self, min_lod_pixels=100, max_lod_pixels=-1, min_fade_extent=200, max_fade_extent=800):
        """
        Builds a kml region.
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

    def build_polygon(self, parent_folder, polygon_style, visibility=True):
        """
        Builds a kml polygon into the parent folder. Polygon is build from netcdf flobal attribute geospatial_bounds.
        :param parent_folder:
        :param polygon_style:
        :return: the input parent folder now containig the polygon and the polygon itself if that is desired instead.
        """
        logger.debug("Building polygon...")
        # get the polygon points from the netcdf global attributes
        try:
            polygon_bounds = self.polygon
            logger.debug(polygon_bounds)
            polygon_bounds = re.sub('POLYGON\(\(', '', polygon_bounds)  # cut out polygon and opening brackets
            polygon_bounds = re.sub('\)\)', '', polygon_bounds)  # cut out trailing brackets
            polygon_bounds = polygon_bounds.split(',')  # turn the string into a list, seperating on the commas
            polygon_bounds = [tuple(p.split(' ')) for p in
                              polygon_bounds]  # within each coord set split the lat and long - group the set in a tuple

            # build the polygon based on the bounds. Also set the polygon name. It is inserted into the parent_folder.
            pol = parent_folder.newpolygon(name=str(self.survey_title) + " " + str(self.survey_id), outerboundaryis=polygon_bounds, visibility=visibility)


            thredds_metadata_link = "http://dapds00.nci.org.au/thredds/catalog/uc0/rr2_dev/axi547/ground_gravity/point_datasets/catalog.html?dataset=uc0/rr2_dev/axi547/ground_gravity/point_datasets/"
            description_string = '<![CDATA['
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', str(self.survey_title))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', str(self.survey_id))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('NCI Data Link', str(thredds_metadata_link) + str(self.survey_id) + '.nc')
            description_string = description_string + ']]>'

            pol.description = description_string

            pol.style = polygon_style
        except:
            logger.debug("No polygon data found.")

        return parent_folder


    def build_points(self, points_folder, bounding_box, point_style):

        logger.debug("Building points for netcdf file: " + str(self.netcdf_path))

        #self.netcdf_dataset = netCDF4.Dataset(self.netcdf_path)
        t1 = time.time()  # Create NetCDFPointUtils object for specified netCDF dataset

        logger.debug('bounding_box')
        logger.debug(bounding_box)

        t2 = time.time()  # create the spatial mask
        bounding_box_floats = [float(coord) for coord in bounding_box]
        spatial_mask = self.npu.get_spatial_mask(bounding_box_floats)
        logger.debug(spatial_mask)
        t3 = time.time()  # get the points and variable info from point generator
        if True in spatial_mask:
            # when ordered through the all_point_data_generator it appears to come out as [obsno, lat, long, then everything else as in field_list]
            field_list = ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype',
                          'reliab', 'gridflag']  # , 'freeair', '', '']

            new_survey_folder = points_folder.newfolder(name=str(self.survey_title) + " " +
                                                                  str(self.survey_id))

            point_data_generator = self.npu.all_point_data_generator(field_list, spatial_mask)
            logger.debug(point_data_generator)
            variable_attributes = next(point_data_generator)
            logger.debug(variable_attributes)
            skip_points = 1
            points_read = 0
            t4 = time.time()  # loop through the points and create them.
            for point_data in point_data_generator:
                points_read += 1

                # ignore points between skip_points
                if points_read % skip_points != 0:
                    continue

                # add new points with netcdf file Obsno as title and long and lat as coordinatess
                new_point = new_survey_folder.newpoint(name="Observation no. " + str(point_data[0]),
                                                   coords=[(point_data[2], point_data[1])])

                description_string = self.build_html_description_string(field_list, variable_attributes, point_data)
                logger.debug(description_string)
                new_point.description = description_string  # set description to point

                if point_data[8] == "Station not used in the production of GA grids.":  # if gridflag equals 2
                    new_point.style.iconstyle.color = 'ff000000'
                    new_point.style.iconstyle.scale = 0.7
                    new_point.style.labelstyle.scale = 0  # removes the label
                    new_point.style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"

                else:
                    new_point.style.iconstyle.color = self.value2colourhex(point_data[5])
                    new_point.style.iconstyle.scale = 0.7
                    new_point.style.labelstyle.scale = 0  # removes the label
                    new_point.style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"

            t5 = time.time()

            return new_survey_folder
        else:
            #logger.debug("no points in view")
            return None


    def build_html_description_string(self, field_list, variable_attributes, point_data):
        """
            Helper function to build the description string automatically based on how many fields are in the field list.
            It will take into account if a unit of measure needs to be specificed or not.
            :param field_list: The fields to be included in the description.
            :param variable_attributes: a list of lists? Matching point_data.
            :param point_data: A list of the actual data of each point
            :return: return the description string html ready to be attached to a kml point.
        """
        description_string = '<![CDATA['
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', self.survey_title)
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', self.survey_id)

        i = 3  # skip obsno, lat, long in field_list
        while i < len(field_list):
            if variable_attributes[field_list[i]].get('units'):
                description_string = description_string + '<p><b>{0}: </b>{1} {2}</p>'.format(
                    variable_attributes[field_list[i]].get('long_name'),
                    point_data[i],
                    variable_attributes[field_list[i]].get('units'))
            else:
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format(
                    variable_attributes[field_list[i]].get('long_name'),
                    point_data[i])
            i += 1
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
        logger.debug("Density: {}".format(((self.east_extent - self.west_extent) * (self.north_extent - self.south_extent)) / self.npu.point_count * 1000))


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

        #spatial_mask = self.npu.get_spatial_mask([130, -40, 150, -5])

        polygon_folder = self.kml.newfolder(name="polygon")
        polygon_folder, polygon = self.build_polygon(polygon_folder, polygon_style)
        dataset_polygon_region = self.build_region(-1, 600, 200, 800)
        dataset_points_region = self.build_region(0, -1, 200, 800)

        points_folder = self.kml.newfolder(name="points")
        points_folder = self.build_points(points_folder, bbox, point_style)

        # structure them correctly
        polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
        points_folder.region = dataset_points_region  # insert built point region into point folder
        return self.kml


    def value2colourhex(self, data_value):
        '''
        Function to convert data value to hex string for color in self.colormap
        '''
        cmap_value = int(convert_value_from_old_to_new_range(data_value, *COLOUR_STRETCH_RANGE, 0, 255))
        #print(cmap_value)
        #print(str(self.colormap(cmap_value)[:3]))
        #print(str(mpl.colors.rgb2hex(self.colormap(cmap_value)[:3])))
        #print(mpl.colors.rgb2hex(self.colormap(cmap_value)[:3]))
        
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
    kml = simplekml.Kml()
    new_folder = kml.newfolder()
    build_dynamic_network_link(new_folder)
    kml.save("dynamic_grav_surveys.kml")

    # --------------------------------
    # build static kml
    # --------------------------------
    # sdmc = SQLiteDatasetMetadataCache(debug=False)
    # endpoint_list = sdmc.search_dataset_distributions(
    #     keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
    #     protocol='opendap',
    #     ll_ur_coords=[[100, -50], [159, -5]]
    #     )
    #
    # for survey in endpoint_list:
    #     netcdf2kml_obj = NetCDF2kmlConverter(survey)
    #     if netcdf2kml_obj.npu.point_count > 2:
    #         static_kml = netcdf2kml_obj.build_static_kml()
    #         static_kml.save("C:\\Users\\u62231\\Desktop\\grav_kmls\\" + netcdf2kml_obj.survey_title + " " + netcdf2kml_obj.survey_id + ".kml")






    # structure them correctly
    # dataset_points_folder.region = dataset_points_region  # insert built point region into point folder
    #dataset_network_link.region = points_region
    # find_subset()
    #network_link = NetCDF2kmlConverter.build_dynamic_network_link(points_folder)


    # build the polygon, points, network link, and regions
    #build_dynamic_kml()
    # server_url_for_dynamic_kml_generation = sys.argv[1]
    # netcdf_path = sys.argv[2]
    # netcdf_path_2 = sys.argv[2]
    # print(server_url_for_dynamic_kml_generation)
    # print(netcdf_path)
    # print(sys.argv)

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

    #modified_server_url_for_dynamic_kml_generation = args.server + 'query'
    #NetCDF2kmlConverter(args.netcdf_path_list)

if __name__ == '__main__':
    main()
