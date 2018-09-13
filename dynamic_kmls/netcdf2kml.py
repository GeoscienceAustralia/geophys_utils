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
import matplotlib.cm as mpl_cm
import matplotlib.colors as mpl_colors
import logging
from datetime import date, timedelta
from geophys_utils.dataset_metadata_cache import SQLiteDatasetMetadataCache
import netCDF4
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils
import numpy as np
import os
from dynamic_kmls import DEBUG

# Setup logging handlers if required
logger = logging.getLogger(__name__)  # Get logger
logger.setLevel(logging.INFO)  # Initial logging level for this module

# Default values when not specified in settings
DEFAULT_POLYGON_COLOUR = 'B30000ff'
DEFAULT_POLYGON_OUTLINE = 1
DEFAULT_LINE_COLOUR = simplekml.Color.yellowgreen
DEFAULT_LINE_WIDTH = 5
DEFAULT_POINT_ICON_SCALE = 0.7
DEFAULT_POINT_LABELSTYLE_SCALE = 0 # 0 removes the label
DEFAULT_POINT_ICON_HREF = 'http://maps.google.com/mapfiles/kml/paddle/grn-blank.png'
DEFAULT_FLAGGED_POINT_ICON_COLOUR = 'ff000000'
DEFAULT_POINT_COLOUR = 'ff32cd9a'

class NetCDF2kmlConverter(object):
    '''
    NetCDF2kmlConverter class definition
    '''
    def __init__(self, netcdf_path, settings, dataset_type, metadata_dict=None):
        '''
        Constructor for NetCDF2kmlConverter class
        @param netcdf_path: netCDF file path or OPeNDAP endpoint
        @param dataset_settings: Dataset settings as read from netcdf2kml_settings.yml settings file
        @param metadata_dict: Dict containing dataset metadata as returned by DatasetMetadataCache.search_dataset_distributions function
        '''
        logger.debug(metadata_dict)
        self.point_utils = None
        self.survey_id = metadata_dict['ga_survey_id']
        self.survey_title = metadata_dict['dataset_title']
        #self.netcdf_path = metadata_dict['distribution_url']  # test
        self.polygon = metadata_dict['convex_hull_polygon']

        self.kml = simplekml.Kml()
        # Store dataset spatial extents as python variables
        self.west_extent = metadata_dict['longitude_min']
        self.east_extent = metadata_dict['longitude_max']
        self.south_extent = metadata_dict['latitude_min']
        self.north_extent = metadata_dict['latitude_max']

        self.point_count = metadata_dict['point_count']

        self.start_date = metadata_dict['start_date']
        self.end_date = metadata_dict['end_date']
        # Set self.end_date if unknown and self.start_date is known
        if self.start_date and not self.end_date:
            self.end_date = self.start_date + timedelta(days=30)

        self.netcdf_path = str(netcdf_path).strip()

        # Don't connect to netCDF until we have to
        self.netcdf_dataset = None
        self.point_utils = None
        self.line_utils = None
        #self.grid_utils = None
        
        # Global settings
        self.line_segments_across_bbox = settings['line_segments_across_bbox']
        
        
        dataset_settings = settings['dataset_settings'][dataset_type]
        
        self.colormap = mpl_cm.get_cmap(dataset_settings['colormap_name'], 
                                        dataset_settings['color_count'])

        self.dataset_link = dataset_settings.get('dataset_link')
        if self.dataset_link:
            # Perform any substitutions required
            self.dataset_link = self.dataset_link.replace('{nc_basename}', os.path.basename(netcdf_path))
            for key, value in metadata_dict.items():
                self.dataset_link = self.dataset_link.replace('{'+key+'}', str(value))
        
        self.polygon_style = simplekml.Style()
        self.polygon_style.polystyle.color = dataset_settings.get('polygon_color') or DEFAULT_POLYGON_COLOUR
        self.polygon_style.polystyle.outline = dataset_settings.get('polygon_outline') or DEFAULT_POLYGON_OUTLINE
        
        self.line_colour = dataset_settings.get('line_colour') or DEFAULT_LINE_COLOUR
        self.line_width = dataset_settings.get('line_width') or DEFAULT_LINE_WIDTH

        self.point_filter = dataset_settings.get('point_filter') # Dict containing field names & values to filter out

        self.point_style = simplekml.Style()
        self.point_style.iconstyle.scale = dataset_settings.get('point_icon_scale') or DEFAULT_POINT_ICON_SCALE
        self.point_style.labelstyle.scale = dataset_settings.get('point_labelstyle_scale') or DEFAULT_POINT_LABELSTYLE_SCALE
        self.point_style.iconstyle.icon.href = dataset_settings.get('point_icon_href') or DEFAULT_POINT_ICON_HREF
        
        self.filtered_point_icon_color = DEFAULT_FLAGGED_POINT_ICON_COLOUR
        
        point_colour_settings = dataset_settings.get('point_colour')
        if point_colour_settings:
            if type(point_colour_settings) == list: # Variable point colour based on specified field value and range
                self.point_colour_field = point_colour_settings[0]
                self.point_colour_range = point_colour_settings[1]
                self.point_colour = None
            elif type(point_colour_settings) == str: # Fixed colour defined in settings
                self.point_colour_field = None
                self.point_colour_range = None
                self.point_colour = point_colour_settings
            else:
                raise BaseException('Invalid point colour settings')
        else: # No colour mapping defined - use default colour
            self.point_colour_field = None
            self.point_colour_range = None
            self.point_colour = DEFAULT_POINT_COLOUR
        
        self.field_list = dataset_settings.get('point_field_list')
        
        self.height_variable = dataset_settings.get('height_variable')


    def __del__(self):
        '''
        NetCDF2kmlConverter destructor
        '''
        if self.netcdf_dataset:
            self.netcdf_dataset.close()
    
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
        # get the polygon points from the polygon string
        try:
            if self.polygon:
                polygon_bounds = [[float(ordinate)
                                   for ordinate in coord_pair.strip().split(' ')
                                   ]
                                  for coord_pair in
                                  re.search('POLYGON\(\((.*)\)\)',
                                            self.polygon
                                            ).group(1).split(',')
                                  ]
                # build the polygon based on the bounds. Also set the polygon name. It is inserted into the parent_folder.
                dataset_folder = parent_folder.newpolygon(name=str(self.survey_title) + " " + str(self.survey_id),
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
                if self.dataset_link:
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Data Link', str(
                    self.dataset_link))
                description_string = description_string + ']]>'
                dataset_folder.description = description_string
    
                dataset_folder.style = self.polygon_style
                
                self.set_timestamps(dataset_folder)
    
                return dataset_folder
        
        except Exception as e:
            #logger.warning("Unable to display polygon "{}": {}".format(self.polygon, e))
            pass

    

    def build_lines(self, parent_folder, bounding_box):
        
        self.netcdf_dataset = self.netcdf_dataset or netCDF4.Dataset(self.netcdf_path)
        self.line_utils = self.line_utils or NetCDFLineUtils(self.netcdf_dataset, 
                                                             enable_disk_cache=False,
                                                             enable_memory_cache=True,
                                                             debug=DEBUG)
        self.point_utils = self.line_utils # NetCDFLineUtils is a subclass of NetCDFPointUtils
        
        #=======================================================================
        # # Set up colour map for rainbow colour scheme
        # line_index_range = (0, len(self.line_utils.line)-1)
        # line_colour_map = {self.line_utils.line[line_index]: self.value2colourhex(line_index, 
        #                                                                      line_index_range)
        #                    for line_index in range(0, len(self.line_utils.line))
        #                    }
        #=======================================================================

        logger.debug("Building lines...")
        bounding_box_floats = [float(coord) for coord in bounding_box]
        
        # Compute segment length as a proportion of the height of bounding box
        subsampling_distance = (bounding_box_floats[3] - bounding_box_floats[1]) / self.line_segments_across_bbox
        
        if self.height_variable:
            height_variable = self.height_variable # e.g. 'lidar'
        else:
            height_variable = [] # Empty string to return no variables, just 'coordinates'
        
        dataset_folder = parent_folder.newfolder(name=str(self.survey_title) + " " + str(self.survey_id))
        
        for line_number, line_data in self.line_utils.get_lines(line_numbers=None, 
                                                                variables=height_variable, 
                                                                bounds=bounding_box_floats,
                                                                subsampling_distance=subsampling_distance
                                                                ):
            #logger.debug("line_number: {}".format(line_number))
            #logger.debug("line_data: {}".format(line_data))
            points_in_subset = len(line_data['coordinates'])
            if points_in_subset:
                line_string = dataset_folder.newlinestring(name=str("Line number: {}".format(line_number)))

                if self.height_variable: # 3D
                    subset_array = np.zeros(shape=(points_in_subset, 3), dtype=line_data['coordinates'].dtype)
                    # Populate coords_3d_array with (x,y,z) coordinates
                    subset_array[:,0:2] = line_data['coordinates']      
                    subset_array[:,2] = line_data[height_variable] # Height above ground
                    
                    line_string.altitudemode = simplekml.AltitudeMode.relativetoground
                else: # 2D
                    subset_array = line_data['coordinates']
                    line_string.altitudemode = simplekml.AltitudeMode.clamptoground
                    
                line_string.coords = subset_array

                line_string.extrude = 0
                line_string.tessellate = 1
                line_string.style.linestyle.width = self.line_width
                
                line_string.style.linestyle.color = self.line_colour # Fixed colour
                #line_segment.style.linestyle.color = line_colour_map[line_number] # Rainbow colour scheme for lines
                
                # build the line description
                description_string = '<![CDATA['
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ', str(self.survey_title))
                description_string = description_string + ']]>'
                line_string.description = description_string

                # Uncomment the following line to enable timestamps on lines
                #self.set_timestamps(line_segment)

               
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
                logger.debug("line doesn't have any points in view")

        return dataset_folder

    def build_points(self, parent_folder, bounding_box):
        """
        Builds all points for a survey. Including building the containing folder, setting time stamps, setting the
         style, and setting the description html to pop up when the point is selected.
        :param parent_folder: The folder for the new survey folder containing the points to be inserted into.
        :return: the kml folder of the survey containing the new points.
        """
        logger.debug("Building points for netcdf file: " + str(self.netcdf_path))
        logger.debug('bounding_box:' + str(bounding_box))
        
        self.netcdf_dataset = self.netcdf_dataset or netCDF4.Dataset(self.netcdf_path)
        self.point_utils = self.point_utils or NetCDFPointUtils(self.netcdf_dataset, 
                                                                enable_disk_cache=False, 
                                                                enable_memory_cache=True,
                                                                debug=DEBUG)
        
        if not self.point_utils.point_count: # No points in dataset
            return None

        bounding_box_floats = [float(coord) for coord in bounding_box]
        logger.debug(bounding_box_floats)
        # bounding_box_floats = [110.4899599829594, -56.11642075733719, 166.658146968822, -6.11642075733719]
        spatial_mask = self.point_utils.get_spatial_mask(bounding_box_floats)
        logger.debug(spatial_mask)
        if np.any(spatial_mask):
            logger.debug("TRUE")
            dataset_folder = parent_folder.newfolder(name=str(self.survey_title) + " " + str(self.survey_id))

            # Set timestamp
            # start_date = re.match('^[0-9]{4}', str(self.survey_id)).group()
            # dataset_folder.timespan.begin = str(start_date) + "-01-01"
            # dataset_folder.timespan.end = str(start_date) + "-01-01"

            point_data_generator = self.point_utils.all_point_data_generator(self.field_list, spatial_mask)
            logger.debug(point_data_generator)
            variable_attributes = next(point_data_generator)
            logger.debug(variable_attributes)
            logger.debug("variable_attributes: " + str(variable_attributes))

            skip_points = 1  # set to limit the points displayed if required.
            points_read = 0

            for point_data_list in point_data_generator:
                point_data = dict(zip(self.field_list, point_data_list)) # Create dict for better readability
                logger.debug("POINT DATA: {}".format(point_data))
                points_read += 1

                # ignore points between skip_points
                if points_read % skip_points != 0:
                    continue

                # add new points with netcdf file Obsno as title and long and lat as coordinatess
                # point_field_list: ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype', 'reliab', 'gridflag']
                new_point = dataset_folder.newpoint(name="Observation no. " + str(point_data['obsno']),
                                                       coords=[(point_data['longitude'], point_data['latitude'])])

                new_point.style = simplekml.Style()
                new_point.style.iconstyle.scale = self.point_style.iconstyle.scale
                new_point.style.labelstyle.scale = self.point_style.labelstyle.scale
                new_point.style.iconstyle.icon.href = self.point_style.iconstyle.icon.href
                
                description_string = self.build_html_description_string(variable_attributes, point_data)
                logger.debug(description_string)
                new_point.description = description_string  # set description to point
                
                # Uncomment the following line to enable timestamps on points
                #self.set_timestamps(new_point)

                # styling
                if self.point_colour: # Fixed point colour
                    unfiltered_point_icon_colour = self.point_colour
                else: # Variable point colour
                    unfiltered_point_icon_colour = self.value2colourhex(point_data[self.point_colour_field], self.point_colour_range)
                
                new_point.style.iconstyle.color = unfiltered_point_icon_colour

                # Set the style for filtered points
                if self.point_filter:  # if there is a point_flag separate the points and colour differently
                    logger.debug('self.point_filter: {}'.format(self.point_filter))
                    for key, value in self.point_filter.items():
                        if point_data[key] == value:
                            new_point.style.iconstyle.color = self.filtered_point_icon_color

            dataset_folder.region = self.build_region(100, -1, 200, 800)
            return dataset_folder


    def build_html_description_string(self, variable_attributes, point_data):
        """
            Helper function to build the description string automatically based on how many fields are in the field list.
            It will take into account if a unit of measure needs to be specified or not.
            :param field_list: The fields to be included in the description.
            :param variable_attributes: a list of lists? Matching self.field_list.
            :param point_data: A dict of the actual data of each point
            :return: return the description string html ready to be attached to a kml point.
        """

        logger.debug(variable_attributes)

        description_string = '<![CDATA['
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', self.survey_title)
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', self.survey_id)

        logger.debug("self.field_list: {}".format(self.field_list))
        for field in self.field_list:
            # skip obsno, lat, long in field_list
            if field in ['obsno', 'latitude', 'longitude']:
                continue
            
            logger.debug(field)
            logger.debug(point_data[field])
            if variable_attributes[field].get('units'):
                description_string = description_string + '<p><b>{0}: </b>{1} {2}</p>'.format(
                    variable_attributes[field].get('long_name'),
                    point_data[field],
                    variable_attributes[field].get('units'))
            else:
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format(
                    variable_attributes[field].get('long_name'),
                    point_data[field])

        description_string = description_string + ']]>'
        return description_string

    #===========================================================================
    # def get_basic_density_measurement(self):
    #     """
    #     Basic calculation to get an idea of the dataset point density. This could be used to dynamically set the
    #     min_lod_pixels, max_lod_pixels, and point scale to get the best visual outcome. A higher density may need a smaller
    #     point scale or even an additional region.
    #     :return: Currently just prints to output
    #     """
    #     logger.debug("Point Count: {}".format(self.point_utils.point_count))
    #     logger.debug("Area: {}".format((self.east_extent - self.west_extent) * (self.north_extent - self.south_extent)))
    #     logger.debug("Density: {}".format(((self.east_extent - self.west_extent) * (
    #     self.north_extent - self.south_extent)) / self.point_utils.point_count * 1000))
    #===========================================================================

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
        dataset_folder = kml.newfolder(name=self.survey_title + " " + self.survey_id)

        dataset_folder = self.kml.newfolder(name="polygon")
        dataset_folder, polygon = self.build_polygon(dataset_folder)
        dataset_polygon_region = self.build_region(-1, 600, 200, 800)
        dataset_points_region = self.build_region(0, -1, 200, 800)

        points_folder = dataset_folder.newfolder(name="points")
        points_folder = self.build_points(points_folder, bbox)

        # structure them correctly
        dataset_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
        points_folder.region = dataset_points_region  # insert built point region into point folder
        return self.kml


    def value2colourhex(self, data_value, data_range, colormap=None):
        '''
        Function to convert data value to hex string for color in colormap_name
        @param data_value: Data value for which to compute colour
        @param data_range: Tuple defining (min, max) of data values over which to compute colours
        @param colormap: matplotlib.mpl_cm.colormap obnect, or None to use class default colormap
        @param clip: Boolean parameter indicating whether out-of-range values should be clipped to range
        
        @return colourhex: KML colour definition string of form ff<rr><gg><bb>, e.g: 'ff60fac5'
        '''
        colormap = colormap or self.colormap
        
        normalised_value = (data_value - data_range[0]) / (data_range[1] - data_range[0])
        
        return mpl_colors.rgb2hex(colormap(normalised_value)[:3]).replace('#', 'ff')


    def set_timestamps(self, kml_entity):
        '''
        Function to set timestamps on specified kml_entity
        '''
        # Add timestamp
        # if self.start_date is not None and self.end_date is not None:
        try:
            assert self.start_date > date(1900, 1, 1), 'Start date {} is less than 1900-01-01'.format(
                self.start_date)
            assert self.end_date < date(2020, 1, 1), 'End date {} is greater than 2020-01-01'.format(self.end_date)
            assert self.end_date > date(1900, 1, 1), 'End date {} is less than 1900-01-01'.format(self.end_date)
            assert self.start_date < date(2020, 1, 1), 'Start date {} is greater than 2020-01-01'.format(
                self.start_date)
            kml_entity.timespan.begin = str(self.start_date)
            kml_entity.timespan.end = str(self.end_date)

        except:  # if survey does not contain start/end date information, use survey id as start/end date year.
            if self.survey_id:
                survey_year = int(re.match('^[0-9]{4}', str(self.survey_id)).group())
                assert survey_year > 1900 and survey_year < 2020, 'survey_year <= 1900 or survey_year >= 2020'

                kml_entity.timespan.begin = str(survey_year) + "-06-01"
                kml_entity.timespan.end = str(survey_year) + "-07-01"


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


#===============================================================================
# def main():
#     # --------------------------------
#     #   build dynamic_grav kml
#     # --------------------------------
#     # kml = simplekml.Kml()
#     # new_folder = kml.newfolder()
#     # build_dynamic_network_link(new_folder)
#     # kml.save("dynamic_grav_surveys.kml")
# 
#     # --------------------------------
#     # build static kml
#     # --------------------------------
#     bounding_box = [[100.00, -50.00], [159.00, -5.00]]  # include all of aus
#     sdmc = SQLiteDatasetMetadataCache(debug=False)
#     endpoint_list = sdmc.search_dataset_distributions(
#         keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
#         protocol='opendap',
#         ll_ur_coords=bounding_box
#     )
#     for survey in endpoint_list:
#         kml = simplekml.Kml()
#         dataset_folder = kml.newfolder(name="Ground Gravity Survey Observations")
#         netcdf2kml_obj = NetCDF2kmlConverter(survey)
#         netcdf2kml_obj.netcdf_dataset = netCDF4.Dataset(netcdf2kml_obj.netcdf_path)
#         netcdf2kml_obj.npu = NetCDFPointUtils(netcdf2kml_obj.netcdf_dataset)
# 
#         survey_region = netcdf2kml_obj.build_region()
#         polygon_region = netcdf2kml_obj.build_region(0, 100)
# 
#         dataset_folder = kml.newfolder(name="polygon")
#         dataset_folder = netcdf2kml_obj.build_polygon(dataset_folder)
#         dataset_folder.region = polygon_region  # insert built point region into point folder
# 
#         dataset_folder.region = survey_region  # insert built point region into point folder
# 
#         if netcdf2kml_obj.npu.point_count > 2:
# 
#             survey_points_folder = netcdf2kml_obj.build_points(dataset_folder,
#                                                                ['110.00', '-45.00', '155.00', '-10.00'])
#             survey_points_folder.region = survey_region  # insert built point region into point folder
# 
#             cleaned_survey_title = re.sub("/", " or ", netcdf2kml_obj.survey_title)
#             cleaned_survey_title = re.sub('"', "", cleaned_survey_title)
# 
#             kml.save('C:\\Users\\u62231\\Desktop\\grav_kmls\\' + str(cleaned_survey_title) + ' ' +
#                      str(netcdf2kml_obj.survey_id) + '.kml')
#         else:
#             logger.debug("fail")
#             netcdf2kml_obj.netcdf_dataset.close()
# 
# 
# 
#             # netcdf_path_list = []
#             # i = 2
#             # while i < len(sys.argv):
#             #     logger.debug(sys.argv[i])
#             #     netcdf_path_list.append(sys.argv[i])
#             #     i += 1
#             # logger.debug(netcdf_path_list)
# 
#             # parser = argparse.ArgumentParser()
#             #
#             # parser.add_argument("-s", "--server", help="The server to receive the get request from google earth and dynamically "
#             #                                            "build kml points within the bbox. If this parameter is empty, a static "
#             #                                            "kml will be generated", type=str, required=False)
#             # parser.add_argument("-n", "--netcdf_path_list", help="Add one or more paths to netcdf files to be converted into a"
#             #                                                      "single kml file.", type=str, nargs='+')
#             #
#             # args = parser.parse_args()
#             #
#             # logger.debug(args.server)
#             # logger.debug(args.netcdf_path_list)
# 
#             # if args.server:
#             #     # dynamic build
#             #     logger.debug("dynamic build")
#             #     if len(args.netcdf_path_list) > 1:
#             #         logger.debug("multiple surveys")
#             #         # multiples
#             #         list_of_converter_objects= []
#             #         for netcdf in args.netcdf_path_list:
#             #             #list_of_converter_objects.append(NetCDF2kmlConverter(netcdf))
#             #             pass
#             #
#             #         # then add the network link using this args.server
# 
#             #     else:
#             #         # single
#             #         logger.debug("single survey")
#             #         converter_object = NetCDF2kmlConverter(args.netcdf_path_list[0])
#             #         converter_object.build_dynamic_kml()
#             #         converter_object.kml.save(converter_object.survey_title + " dynamic points.kml")
#             #         logger.debug("Building kml for survey: " + converter_object.survey_title + " dynamic points.kml")
#             # else:
#             #     # static build
#             #     logger.debug("static build")
#             #     if len(args.netcdf_path_list) > 1:
#             #         logger.debug("multiple surveys")
#             #         pass
#             #     else:
#             #         logger.debug("single survey")
#             #         converter_object = NetCDF2kmlConverter(args.netcdf_path_list[0])
#             #         converter_object.build_static_kml()
#             #
#             #         converter_object.kml.save(converter_object.survey_title + " static points.kml")
#             #         logger.debug("Building kml for survey: " + converter_object.survey_title + " static points.kml")
#             # single
# 
#             # modified_server_url_for_dynamic_kml_generation = args.server + 'query'
#             # NetCDF2kmlConverter(args.netcdf_path_list)
# 
# 
# if __name__ == '__main__':
#     main()
#===============================================================================
