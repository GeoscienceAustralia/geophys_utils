import simplekml
import netCDF4
import numpy as np
from geophys_utils import NetCDFPointUtils
import re
import sys
import argparse
import time


class NetCDF2kmlConverter(object):

    def __init__(self, netcdf_path):
        self.netcdf_path = netcdf_path

        # Create NetCDFPointUtils object for specified netCDF dataset
        #netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc'
        #netcdf_path = 'E:\\Temp\\gravity_point_test\\195256.nc'
        #195256
        #netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_netcdf_4\\201780.nc"
        #195105
        #201780

        self.netcdf_dataset = netCDF4.Dataset(self.netcdf_path)
        self.npu = NetCDFPointUtils(self.netcdf_dataset)
        self.survey_id = str(self.npu.netcdf_dataset.getncattr('survey_id'))
        self.survey_title = str(self.npu.netcdf_dataset.getncattr('title'))
        self.kml = simplekml.Kml()

        # Store dataset spatial extents as python variables
        self.west_extent = self.npu.netcdf_dataset.getncattr('geospatial_lon_min')
        self.east_extent = self.npu.netcdf_dataset.getncattr('geospatial_lon_max')
        self.south_extent = self.npu.netcdf_dataset.getncattr('geospatial_lat_min')
        self.north_extent = self.npu.netcdf_dataset.getncattr('geospatial_lat_max')
        print("WEST EXTENT")
        print(self.west_extent)
        # set kml region constants
        # Measurement in screen pixels that represents the maximum limit of the visibility range for a given Region.
        self.MIN_LOD_PIXELS = 1200
        self.MAX_LOD_PIXELS = -1  # -1 the default, indicates "active to infinite size."

        # Distance over which the geometry fades, from fully opaque to fully transparent.
        # This ramp value, expressed in screen pixels, is applied at the minimum end of the LOD (visibility) limits.


        # point constants
        self.POINT_ICON_STYLE_LINK = "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"

        # print(west_extent)
        # print("MASK")
        # mask = np.zeros(shape=(npu.point_count,), dtype='bool')
        # mask[-10:] = True
        # print(mask)

    # -----------------------bounds------------
    #        self.bounds = [xmin, ymin, xmax, ymax]

    def get_basic_density_measurement(self):
        """
        Basic calculation to get an idea of the dataset point density. This could be used to dynamically set the
        min_lod_pixels, max_lod_pixels, and point scale to get the best visual outcome. A higher density may need a smaller
        point scale or even an additional region.
        :return: Currently just prints to output
        """
        print("Point Count: {}".format(self.npu.point_count))
        print("Area: {}".format((self.east_extent - self.west_extent) * (self.north_extent - self.south_extent)))
        print("Density: {}".format(((self.east_extent - self.west_extent) * (self.north_extent - self.south_extent)) / self.npu.point_count * 1000))


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


    def build_polygon(self, parent_folder, polygon_style):
        """
        Builds a kml polygon into the parent folder. Polygon is build from netcdf flobal attribute geospatial_bounds.
        :param parent_folder:
        :param polygon_style:
        :return: the input parent folder now containig the polygon and the polygon itself if that is desired instead.
        """

        # get the polygon points from the netcdf global attributes
        polygon_bounds = re.sub(',\s', ',', self.npu.netcdf_dataset.geospatial_bounds) # cut out space after comma between coord sets
        polygon_bounds = re.sub('POLYGON\(\(', '', polygon_bounds) # cut out polygon and opening brackets
        polygon_bounds = re.sub('\)\)', '', polygon_bounds) # cut out trailing brackets
        polygon_bounds = polygon_bounds.split(',') # turn the string into a list, seperating on the commas
        polygon_bounds = [tuple(p.split(' ')) for p in polygon_bounds] # within each coord set split the lat and long - group the set in a tuple

        # build the polygon based on the bounds. Also set the polygon name. It is inserted into the parent_folder.
        pol = parent_folder.newpolygon(name=str(self.survey_title) + " " + str(self.survey_id),
                             outerboundaryis=polygon_bounds)

        pol.style = polygon_style

        return parent_folder, pol


    def plot_points_in_kml(self):
        """
        Loop through points in the netcdf file, and build and style each point as a kml place mark point.
        :return: The point folder. This folder is called to insert the corresponding region
        """
        index = 0
        count = self.npu.point_count
        points_folder = self.kml.newfolder()
        while index < count:
            # add new points with netcdf file Obsno as title and long and lat as coordinatess
            new_point = points_folder.newpoint(name=str(self.npu.netcdf_dataset.variables["obsno"][index]),
            coords=[(self.npu.netcdf_dataset.variables["longitude"][index], self.npu.netcdf_dataset.variables["latitude"][index])])

            # set the point icon. Different urls can be found in point style options in google earth
            new_point.style.iconstyle.icon.href = self.POINT_ICON_STYLE_LINK
            new_point.style.iconstyle.scale = 0.7
            new_point.labelstyle.scale = 0  # removes the label

            index = index + 1
        assert index == count
        return points_folder


    def build_dynamic_network_link(self, containing_folder):
        """
        Build a network link, set the parameters, and inserts into the specified containing folder.
        """
        net_link = containing_folder.newnetworklink(name="Network Link")
        net_link.link.href = "http://127.0.0.1:5000/query"
        net_link.link.viewrefreshmode = simplekml.ViewRefreshMode.onstop
        net_link.link.viewrefreshtime = 1
        net_link.link.refreshinterval = 2

        return net_link



       def find_subset(self):

        lats = self.npu.netcdf_dataset.variables['latitude'][:]
        longs = self.npu.netcdf_dataset.variables['longitude'][:]
        #print(longs)
        latselect = np.logical_and(lats > 38.03, lats < 40)
        #print(latselect)
        lonselect = np.logical_and(longs > 146.6, longs < 147)
        #print(lonselect)

        #print([i for i in latselect if i is True])
        data = self.npu.netcdf_dataset.variables['obsno']
        #data = npu.netcdf_dataset.variables['Obsno'][0, 0, latselect, lonselect]
        #print(data)


    def build_static_kml(self):
        polygon_folder = self.kml.newfolder()
        polygon_folder, polygon = self.build_polygon(polygon_folder)
        dataset_polygon_region = self.build_region(-2, 600, 200, 800)
        dataset_points_region = self.build_region(0, -1, 200, 800)

        points_folder = self.build_points_2()

        # structure them correctly
        polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
        points_folder.region = dataset_points_region  # insert built point region into point folder

        return self.kml


    def build_points(self, points_folder, bounding_box, point_style):

        t1 = time.time()  # Create NetCDFPointUtils object for specified netCDF dataset

        #POINT_ICON_STYLE_LINK = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"

        survey_title = str(self.npu.netcdf_dataset.getncattr('title'))
        print('bounding_box')
        print(bounding_box)

        t2 = time.time()  # create the spatial mask
        bounding_box_floats = [float(coord) for coord in bounding_box]
        spatial_mask = self.npu.get_spatial_mask(bounding_box_floats)
        t3 = time.time()  # get the points and variable info from point generator
        if True in spatial_mask:
            # when ordered through the all_point_data_generator it appears to come out as [obsno, lat, long, then everything else as in field_list]
            field_list = ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype',
                          'reliab', 'gridflag']  # , 'freeair', '', '']

            point_data_generator = self.npu.all_point_data_generator(field_list, spatial_mask)
            variable_attributes = next(point_data_generator)

            skip_points = 1
            points_read = 0
            t4 = time.time()  # loop through the points and create them.
            for point_data in point_data_generator:
                points_read += 1

                # ignore points between skip_points
                if points_read % skip_points != 0:
                    continue

                # add new points with netcdf file Obsno as title and long and lat as coordinatess
                new_point = points_folder.newpoint(name="Observation no. " + str(point_data[0]),
                                                   coords=[(point_data[2], point_data[1])])
                # set point html for pop up description.
                description_string = '<![CDATA[' \
                                     '<p><b>{0}: </b>{1} {2}</p>' \
                                     '<p><b>{3}: </b>{4} {5}</p> ' \
                                     '<p><b>{6}: </b>{7} {8}</p>' \
                                     '<p><b>{9}: </b>{10}</p> ' \
                                     '<p><b>{11}: </b>{12}</p>' \
                                     '<p><b>{13}: </b>{14}</p>' \
                                     ']]>'.format(
                    variable_attributes['grav'].get('long_name'), point_data[3],
                    variable_attributes['grav'].get('units'),
                    variable_attributes['freeair'].get('long_name'), point_data[4],
                    variable_attributes['freeair'].get('units'),  # free air
                    variable_attributes['bouguer'].get('long_name'), point_data[5],
                    variable_attributes['bouguer'].get('units'),  # bouguer
                    variable_attributes['stattype'].get('long_name'), point_data[6],  # station type
                    variable_attributes['reliab'].get('long_name'), point_data[7],  # reliability
                    variable_attributes['gridflag'].get('long_name'), point_data[8]  # Gridflag
                )

                new_point.description = description_string  # set description to point
                new_point.style = point_style  # set style to point

            t5 = time.time()


            return points_folder
        else:
            #print("no points in view")
            return None
    # def build_html_description_string(self, field_list, point_data_generator):
    #
    #     for field in field_list:
    #         if variable_attributes

def main():

    # kml = simplekml.Kml()
    # points_folder = kml.newfolder()
    # #dataset_network_link = NetCDF2kmlConverter.build_dynamic_network_link(points_folder)
    #
    # net_link = points_folder.newnetworklink(name="Network Link")
    # net_link.link.href = "http://127.0.0.1:5000/query"
    # net_link.link.viewrefreshmode = simplekml.ViewRefreshMode.onstop
    # net_link.link.viewrefreshtime = 1
    # net_link.link.refreshinterval = 2
    # kml.save("dynamic points.kml")

    netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc'
    static = NetCDF2kmlConverter(netcdf_path)
    static_kml = static.build_static_kml()
    static_kml.save("static_points.kml")





    # structure them correctly
    # dataset_points_folder.region = dataset_points_region  # insert built point region into point folder
    #dataset_network_link.region = points_region
    # find_subset()
    #network_link = NetCDF2kmlConverter.build_dynamic_network_link(points_folder)




    #"Poorly controlled data which should be used cautiously."
    #"Data with weak gravity, position and elevation control."
    ##new_point.description = "description {0} here {1}".format(point_data[2], point_data[3])
    #get_Var = npu.get_lookup_mask(lookup_value_list=[] ,lookup_variable_name='reliab')
    #print(get_Var)
    # build the polygon, points, network link, and regions
    #build_dynamic_kml()
    # server_url_for_dynamic_kml_generation = sys.argv[1]
    # netcdf_path = sys.argv[2]
    # netcdf_path_2 = sys.argv[2]
    # print(server_url_for_dynamic_kml_generation)
    # print(netcdf_path)
    # print(sys.argv)
    #
    #
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
    #
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
    #
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

