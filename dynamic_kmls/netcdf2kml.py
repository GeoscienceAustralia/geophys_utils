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

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
import simplekml
import re
import matplotlib.cm as mpl_cm
import matplotlib.colors as mpl_colors
import logging
from datetime import date, timedelta
import netCDF4
import numpy as np
import os
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils, NetCDFGridUtils

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class NetCDF2kmlConverter(object):
    '''
    NetCDF2kmlConverter class definition
    '''
    def __init__(self, settings, dataset_type, debug=False):
        '''
        Constructor for NetCDF2kmlConverter class
        @param settings: Dataset settings as read from netcdf2kml_settings.yml settings file
        '''
        # Initialise and set debug property
        self._debug = None
        self.debug = debug
        
        logger.debug('Instantiating NetCDF2kmlConverter object for {} datasets'.format(dataset_type))
        
        # Create combined full settings dict and set instance attributes
        combined_settings = dict(settings['global_settings']) # Global settings
        combined_settings.update(settings['default_dataset_settings']) # Default dataset settings
        combined_settings.update(settings['dataset_settings'][dataset_type]) # Dataset settings
        combined_settings['dataset_type'] = dataset_type
        
        logger.debug('combined_settings: {}'.format(combined_settings))        
        for key, value in combined_settings.items():
            setattr(self, key, value)
            
        # Dict of KML buiding functions keyed by dataset form
        self.build_kml_functions = {'polygon': self.build_polygon,
                                    'point': self.build_points,
                                    'line': self.build_lines,
                                    'grid': self.build_thumbnail
                                    }
        
        # Initialise private instance values to None - will be set by property getter methods as required
        self._netcdf_path = None
        self._netcdf_dataset = None
        self._point_utils = None
        self._line_utils = None
        self._grid_utils = None

        self.colormap = mpl_cm.get_cmap(self.colormap_name, 
                                        self.color_count)

        if type(self.point_color_settings) == list: # Variable point color based on specified field value and range
            self.point_color_field = self.point_color_settings[0]
            self.point_color_range = self.point_color_settings[1]
            self.point_color = None
        elif type(self.point_color_settings) == str: # Fixed color defined in settings
            self.point_color_field = None
            self.point_color_range = None
            self.point_color = self.point_color_settings
        else:
            raise BaseException('Invalid point color settings') # Should never happen
                        
        if self.filtered_point_icon_color:
            self.filtered_point_style = simplekml.Style()
            self.filtered_point_style.iconstyle.scale = self.point_icon_scale
            self.filtered_point_style.labelstyle.scale = self.point_labelstyle_scale
            self.filtered_point_style.iconstyle.icon.href = self.point_icon_href
            self.filtered_point_style.iconstyle.color = self.filtered_point_icon_color
        else:
            self.filtered_point_style = None

        self.kml = simplekml.Kml()
        self.dataset_type_folder = self.kml.newfolder(name="No {} in view".format(self.dataset_type_name))
        #=======================================================================
        # # Set fixed styles for sub-elements
        # self.dataset_type_folder.style.polystyle.color = self.polygon_color
        # self.dataset_type_folder.style.polystyle.outline = self.polygon_outline
        # self.dataset_type_folder.style.iconstyle.scale = self.point_icon_scale
        # self.dataset_type_folder.style.labelstyle.scale = self.point_labelstyle_scale
        # self.dataset_type_folder.style.iconstyle.icon.href = self.point_icon_href
        # if self.point_color: # Fixed point colour - will None if variant color required
        #     self.dataset_type_folder.style.iconstyle.color = self.point_color
        # self.dataset_type_folder.style.linestyle.width = self.line_width
        # self.dataset_type_folder.style.linestyle.color = self.line_color # Fixed color
        #=======================================================================
        
        self.polygon_style = simplekml.Style()
        self.polygon_style.polystyle.color = self.polygon_color
        self.polygon_style.polystyle.outline = self.polygon_outline

        self.line_style = simplekml.Style()
        self.line_style.linestyle.width = self.line_width
        self.line_style.linestyle.color = self.line_color # Fixed color
        
        self.point_style = simplekml.Style()
        self.point_style.iconstyle.scale = self.point_icon_scale
        self.point_style.labelstyle.scale = self.point_labelstyle_scale
        self.point_style.iconstyle.icon.href = self.point_icon_href
        if self.point_color: # Fixed point colour - will None if variant color required
            self.point_style.iconstyle.color = self.point_color
        
        # Initialise point style cache for variant colors to cut down the number of style definitions required
        self.point_style_by_color = {}  
        
        self.dataset_count = 0
        

    def __del__(self):
        '''
        NetCDF2kmlConverter destructor
        '''
        if self._netcdf_dataset:
            self._netcdf_dataset.close()
    
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

        region = simplekml.Region(latlonaltbox="<north>" + str(self.latitude_max) + "</north>" +
                                               "<south>" + str(self.latitude_min) + "</south>" +
                                               "<east>" + str(self.longitude_max) + "</east>" +
                                               "<west>" + str(self.longitude_min) + "</west>",
                                  lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
                                      "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
                                      "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
                                      "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")
        return region

    def build_polygon(self, bounding_box, visibility=True):
        """
        Builds a kml polygon into the parent folder. Polygon is built from convex_hull_polygon in database search result.
        @param self.dataset_type_folder: KML folder under which to build geometry for line dataset
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        """
        # get the polygon points from the polygon string
        try:
            if self.convex_hull_polygon:
                polygon_bounds = [[float(ordinate)
                                   for ordinate in coord_pair.strip().split(' ')
                                   ]
                                  for coord_pair in
                                  re.search('POLYGON\(\((.*)\)\)',
                                            self.convex_hull_polygon
                                            ).group(1).split(',')
                                  ]
                # build the polygon based on the bounds. Also set the polygon name. It is inserted into the self.dataset_type_folder.
                dataset_kml = self.dataset_type_folder.newpolygon(name=str(self.dataset_title) + " " + str(self.ga_survey_id),
                                               outerboundaryis=polygon_bounds, visibility=visibility)
                
                dataset_kml.style = self.polygon_style
    
                # Always set timestamps on polygons
                self.set_timestamps(dataset_kml)

                # build the polygon description
                description_string = '<![CDATA['
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name',
                                                                                          str(self.dataset_title))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', str(self.ga_survey_id))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Start Date',
                                                                                          str(self.start_date))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey End Date',
                                                                                          str(self.end_date))
                if self.modified_dataset_link:
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
                    self.modified_dataset_link))
                description_string = description_string + ']]>'
                dataset_kml.description = description_string
    
                return dataset_kml
        
        except Exception as e:
            logger.debug('Unable to display polygon "{}": {}'.format(self.convex_hull_polygon, e))
    

    def build_lines(self, bounding_box, visibility=True):
        '''
        Builds a set of kml linestrings into the parent folder. Linestrings are built from dynamically subsampled points in line.
        @param self.dataset_type_folder: KML folder under which to build geometry for line dataset
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        '''
        bounding_box_floats = [float(coord) for coord in bounding_box]
        
        # Compute segment length as a proportion of the height of bounding box
        subsampling_distance = (bounding_box_floats[3] - bounding_box_floats[1]) / self.line_segments_across_bbox
        
        if self.height_variable:
            height_variable = self.height_variable # e.g. 'lidar'
        else:
            height_variable = [] # Empty string to return no variables, just 'coordinates'
        
        dataset_kml = self.dataset_type_folder.newfolder(name=str(self.dataset_title) + " " + str(self.ga_survey_id), visibility=visibility)
        
        #dataset_kml.style = self.dataset_type_folder.style
        #dataset_kml.style = self.line_style
    
        if self.timestamp_detail_view:
            # Enable timestamps on lines
            self.set_timestamps(dataset_kml)
        
        visible_line_count = 0 # Need to iterate through generator to count lines
        for line_number, line_data in self.line_utils.get_lines(line_numbers=None, 
                                                                variables=height_variable, 
                                                                bounds=bounding_box_floats,
                                                                subsampling_distance=subsampling_distance
                                                                ):
            #logger.debug("line_number: {}".format(line_number))
            #logger.debug("line_data: {}".format(line_data))
            points_in_subset = len(line_data['coordinates'])
            if points_in_subset:
                visible_line_count += 1
                line_string = dataset_kml.newlinestring(name=str("Line number: {}".format(line_number)))
                
                #line_string.style = dataset_kml.style

                line_string.style = self.line_style               
                
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
                
                # build the line description
                description_string = '<![CDATA['
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ', str(self.dataset_title))
                if self.modified_dataset_link:
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
                    self.modified_dataset_link))
                description_string = description_string + ']]>'
                line_string.description = description_string
               
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

        if not visible_line_count:
            dataset_kml.name = '{} has no lines in view'.format(self.dataset_title)
            dataset_kml.visibility = False
            
        return dataset_kml

    def build_points(self, bounding_box, visibility=True):
        """
        Builds all points for a survey. Including building the containing folder, setting time stamps, setting the
         style, and setting the description html to pop up when the point is selected.
        @param self.dataset_type_folder: KML folder under which to build geometry for line dataset
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        """        
        
        visible_points_exist = False
        if self.point_utils.point_count: # Points found in dataset
            spatial_mask = self.point_utils.get_spatial_mask(bounding_box)
            #logger.debug('spatial_mask: {}'.format(spatial_mask))
            visible_points_exist = np.any(spatial_mask)

        if not visible_points_exist:   
            logger.debug('No points in view')
            return
        
        dataset_kml = self.dataset_type_folder.newfolder(name=str(self.dataset_title) + " " + str(self.ga_survey_id), visibility=visibility)
        
        dataset_kml.style = self.point_style

        if self.timestamp_detail_view:
            # Enable timestamps on points
            self.set_timestamps(dataset_kml)
            
        point_data_generator = self.point_utils.all_point_data_generator(self.point_field_list, spatial_mask)
        logger.debug(point_data_generator)
        variable_attributes = next(point_data_generator) # Get point attribute names from first item returned
        logger.debug(variable_attributes)
        logger.debug("variable_attributes: " + str(variable_attributes))

        skip_points = 1  # set to limit the points displayed if required.
        visible_point_count = 0

        for point_data_list in point_data_generator:
            point_data = dict(zip(self.point_field_list, point_data_list)) # Create dict for better readability
            logger.debug("POINT DATA: {}".format(point_data))
            visible_point_count += 1

            # ignore points between skip_points
            if visible_point_count % skip_points != 0:
                continue

            # add new points with netcdf file Obsno as title and long and lat as coordinatess
            # point_field_list: ['obsno', 'latitude', 'longitude', 'grav', 'freeair', 'bouguer', 'stattype', 'reliab', 'gridflag']
            point_kml = dataset_kml.newpoint(name="Observation no. " + str(point_data['obsno']),
                                                   coords=[(point_data['longitude'], point_data['latitude'])])

            point_kml.style = dataset_kml.style
            
            description_string = self.build_html_description_string(variable_attributes, point_data)
            logger.debug(description_string)
            point_kml.description = description_string  # set description of point
            
            # Set point styling               
            # Set the color for filtered points
            variant_point_style = None # Assume point is unfiltered
            if self.point_filter:  # if there is a point_flag separate the points and color differently
                logger.debug('self.point_filter: {}'.format(self.point_filter))
                for key, value in self.point_filter.items():
                    if point_data[key] == value: # Point satisfies filter condition
                        variant_point_style = self.filtered_point_style
                        break
            
            if not variant_point_style: # Point is not filtered 
                if not self.point_color: # Variable point color required
                    variant_point_color = self.value2colorhex(point_data[self.point_color_field], self.point_color_range)
                    
                    variant_point_style = self.point_style_by_color.get(variant_point_color) # Check point style cache
                    if not variant_point_style: # Point style doesn't already exist in cache - construct and cache it
                        variant_point_style = simplekml.Style()
                        variant_point_style.iconstyle.scale = self.point_icon_scale
                        variant_point_style.labelstyle.scale = self.point_labelstyle_scale
                        variant_point_style.iconstyle.icon.href = self.point_icon_href
                        variant_point_style.iconstyle.color = variant_point_color
                        self.point_style_by_color[variant_point_color] = variant_point_style

            # Override default style if required    
            if variant_point_style:
                point_kml.style = variant_point_style

        dataset_kml.region = self.build_region(100, -1, 200, 800)
            
        return dataset_kml


    def build_thumbnail(self, bounding_box, visibility=True):
        """
        Builds a kml thumbnail image into the parent folder. 
        Thumbnail URL is built from OPeNDAP endpoint at this stage, but this needs to change.
        @param self.dataset_type_folder: KML folder under which to build geometry for line dataset
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        """
        logger.debug("Building WMS thumbnail...")
        try:
            logger.debug("Dataset WEST extent: {}".format(self.longitude_min))
            logger.debug("BBOX WEST extent: {}".format(bounding_box[0]))
            logger.debug("Dataset EAST extent: {}".format(self.longitude_max))
            logger.debug("BBOX EAST extent: {}".format(bounding_box[2]))
            logger.debug("Dataset SOUTH extent: {}".format(self.latitude_min))
            logger.debug("BBOX SOUTH extent: {}".format(bounding_box[1]))
            logger.debug("Dataset NORTH extent: {}".format(self.latitude_max))
            logger.debug("BBOX NORTH extent: {}".format(bounding_box[3]))

            # Define smallest bounding box to retrieve via WMS
            west = max(bounding_box[0], self.longitude_min)
            east = min(bounding_box[2], self.longitude_max)
            south = max(bounding_box[1], self.latitude_min)
            north = min(bounding_box[3], self.latitude_max)

            wms_url = self.distribution_url.replace('/dodsC/', '/wms/') #TODO: Replace this hack
            logger.debug("WMS URL")
            logger.debug(self.distribution_url)
            logger.debug(self.metadata_uuid)

            logger.debug(wms_url)
            wms_url = wms_url + "?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&BBOX={0},{1},{2},{3}&CRS=EPSG:4326&WIDTH={4}&HEIGHT={5}&LAYERS={6}&STYLES=&FORMAT=image/png" \
                      "&DPI=120&MAP_RESOLUTION=120&FORMAT_OPTIONS=dpi:120&TRANSPARENT=TRUE" \
                      "&COLORSCALERANGE=-350%2C350&NUMCOLORBANDS=127".format(south, west, north, east, str(int((east - west) / self.wms_pixel_size)), str(int((north - south) / self.wms_pixel_size)), self.wms_layer_name)

            #mag_tmi_anomaly

            # wms_url = "http://dapds00.nci.org.au/thredds/wms/rr2/airborne_geophysics/NSW/P1027/magnetics/grid/mNSW1027/" \
            #           "mNSW1027.nc?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap" \
            #           "&BBOX={0},{1},{2},{3}" \
            #           "&CRS=EPSG:4326&WIDTH=206&HEIGHT=269&LAYERS=mag_tmi_anomaly&STYLES=&FORMAT=image/png" \
            #           "&DPI=120&MAP_RESOLUTION=120&FORMAT_OPTIONS=dpi:120&TRANSPARENT=TRUE" \
            #           "&COLORSCALERANGE=-2781%2C2741&NUMCOLORBANDS=10".format(south, west, north, east)

            # dataset_kml = self.dataset_type_folder.newfolder(name='overlay_test',
            #                                                  visibility=visibility)

            # dataset_kml.style = self.point_style

            dataset_kml = self.dataset_type_folder.newgroundoverlay(name=self.dataset_title)
            logger.debug(wms_url)
            dataset_kml.icon.href = wms_url
            logger.debug(dataset_kml.icon.href)
            # dataset_kml.gxlatlonquad.coords = [(18.410524, -33.903972), (18.411429, -33.904171),
            #                                    (18.411757, -33.902944), (18.410850, -33.902767)]
            logger.debug("NORTH")
            logger.debug(north)
            dataset_kml.latlonbox.north = north
            dataset_kml.latlonbox.south = south
            dataset_kml.latlonbox.east = east
            dataset_kml.latlonbox.west = west
            dataset_kml.color = 'aaffffff'

            logger.debug("GROUND")
            logger.debug(dataset_kml)

            self.set_timestamps(dataset_kml)

            # build the bitmap description
            description_string = '<![CDATA['
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name',
                                                                                      str(self.dataset_title))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', str(self.ga_survey_id))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Start Date',
                                                                                      str(self.start_date))
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey End Date',
                                                                                      str(self.end_date))
            if self.modified_dataset_link:
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
                self.modified_dataset_link))
            description_string = description_string + ']]>'
            
            #TODO: See whether we can make this a displayable balloon from the map
            dataset_kml.description = description_string 

            return dataset_kml
    #===========================================================================
        
        except Exception as e:
            #logger.warning("Unable to display thumbnail "{}": {}".format(wms_url, e))
            pass

    

    def build_html_description_string(self, variable_attributes, point_data):
        """
            Helper function to build the description string automatically based on how many fields are in the field list.
            It will take into account if a unit of measure needs to be specified or not.
            :param field_list: The fields to be included in the description.
            :param variable_attributes: a list of lists? Matching self.point_field_list.
            :param point_data: A dict of the actual data of each point
            :return: return the description string html ready to be attached to a kml point.
        """

        logger.debug(variable_attributes)

        description_string = '<![CDATA['
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', self.dataset_title)
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', self.ga_survey_id)

        logger.debug("self.point_field_list: {}".format(self.point_field_list))
        for field in self.point_field_list:
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

        if self.modified_dataset_link:
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
            self.modified_dataset_link))
            
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
    #     logger.debug("Area: {}".format((self.longitude_max - self.longitude_min) * (self.latitude_max - self.latitude_min)))
    #     logger.debug("Density: {}".format(((self.longitude_max - self.longitude_min) * (
    #     self.latitude_max - self.latitude_min)) / self.point_utils.point_count * 1000))
    #===========================================================================

#===============================================================================
#     def build_static_kml(self):
# 
#         point_style = simplekml.Style()
#         point_style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/paddle/grn-blank.png"
#         point_style.iconstyle.scale = 0.7
#         point_style.labelstyle.scale = 0  # removes the label
# 
#         bbox = ['115.6900403947253', '-80.3741299166767', '165.9763920890655', '-7.21307533348337']
#         # polygon style
#         polygon_style = simplekml.Style()
#         polygon_style.polystyle.color = '990000ff'  # Transparent red
#         polygon_style.polystyle.outline = 1
# 
#         kml = simplekml.Kml()
#         dataset_kml = kml.newfolder(name=self.dataset_title + " " + self.ga_survey_id)
# 
#         dataset_kml = self.kml.newfolder(name="polygon")
#         dataset_kml, polygon = self.build_polygon(dataset_kml)
#         dataset_polygon_region = self.build_region(-1, 600, 200, 800)
#         dataset_points_region = self.build_region(0, -1, 200, 800)
# 
#         points_folder = dataset_kml.newfolder(name="points")
#         points_folder = self.build_points(points_folder, bbox)
# 
#         # structure them correctly
#         dataset_kml.region = dataset_polygon_region  # insert built polygon region into polygon folder
#         points_folder.region = dataset_points_region  # insert built point region into point folder
#         
#         return self.kml
#===============================================================================


    def value2colorhex(self, data_value, data_range, colormap=None):
        '''
        Function to convert data value to hex string for color in colormap_name
        @param data_value: Data value for which to compute color
        @param data_range: Tuple defining (min, max) of data values over which to compute colors
        @param colormap: matplotlib.mpl_cm.colormap obnect, or None to use class default colormap
        @param clip: Boolean parameter indicating whether out-of-range values should be clipped to range
        
        @return colorhex: KML color definition string of form ff<rr><gg><bb>, e.g: 'ff60fac5'
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
            try:
                if self.ga_survey_id:
                    survey_year = int(re.match('^[0-9]{4}', str(self.ga_survey_id)).group())
                    assert survey_year > 1900 and survey_year < 2020, 'survey_year <= 1900 or survey_year >= 2020'
    
                    kml_entity.timespan.begin = str(survey_year) + "-06-01"
                    kml_entity.timespan.end = str(survey_year) + "-07-01"
            except Exception as e:
                logger.debug('Unable to parse year from survey ID: {}'.format(e))
                pass

    def build_dataset_kml(self, kml_format, dataset_metadata_dict, bbox_list, visibility=True):
        '''
        Function to build and return KML of specified format for specified dataset
        @param kml_format: Format of KML required. Must be in ['polygon', 'point', 'line', 'grid']
        @param dataset_metadata_dict: Dict containing dataset metadata as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bbox_list: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        '''
        #logger.debug('build_dataset_kml({}, {}, {}, {}) called'.format(kml_format, dataset_metadata_dict, bbox_list, visibility))
        self.netcdf_path = dataset_metadata_dict['netcdf_path'] # Property setter will take care of dependencies
        
        # Update instance attributes from dataset_metadata_dict
        for key, value in dataset_metadata_dict.items():
            setattr(self, key, value)
            
        #TODO: put this back in self.netcdf_path setter method somehow
        if self.dataset_link:
            # Perform any substitutions required
            dataset_metadata_dict.update({'netcdf_path': self.netcdf_path,
                                          'netcdf_basename': self.netcdf_basename
                                          })
            self.modified_dataset_link = self.dataset_link
            for key, value in dataset_metadata_dict.items():
                self.modified_dataset_link = self.modified_dataset_link.replace('{'+key+'}', str(value))
        else:
            self.modified_dataset_link = None
            
        # Set self.end_date if unknown and self.start_date is known
        if self.start_date and not self.end_date:
            self.end_date = self.start_date + timedelta(days=30)
        
        build_kml_function = self.build_kml_functions.get(kml_format)
        assert build_kml_function, 'Invalid dataset form "{}". Must be in {}'.format(self.dataset_format, 
                                                                                     list(self.build_kml_functions.keys()))
        logger.debug('Processing {}s for dataset {}'.format(self.dataset_format, self.netcdf_path))
        return build_kml_function(bbox_list, visibility)

        
        
    def build_bbox_kml(self, dataset_metadata_dict_list, bbox_list, visibility=True):
        '''
        Function to build and return KML for entire bounding box
        @param dataset_metadata_dict_list: List of dicts containing dataset metadata as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bbox_list: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        '''
        # Show polygons for low zoom
        if (bbox_list[2] - bbox_list[0]) >= self.min_polygon_bbox_width:
            kml_format = 'polygon'
        else:
            kml_format = self.dataset_format
        
        # Reset KML
        self.kml = simplekml.Kml()
        self.dataset_type_folder = self.kml.newfolder(name="No {} in view".format(self.dataset_type_name))

        self.dataset_count = 0
        for dataset_metadata_dict in dataset_metadata_dict_list: 
            # N.B: Could determine visibility from data here
            if self.build_dataset_kml(kml_format, dataset_metadata_dict, bbox_list, visibility):   
                self.dataset_count += 1
        
        if self.dataset_count:
            self.dataset_type_folder.name = '{} {} in view'.format(self.dataset_count, self.dataset_type_name)
            
        
    @property
    def netcdf_dataset(self):
        '''
        Getter method for netcdf_dataset property. Will set self._netcdf_dataset if None
        '''
        if not self._netcdf_dataset:
            # Don't connect to netCDF dataset until we have to
            self._netcdf_dataset = netCDF4.Dataset(self.netcdf_path)
        return self._netcdf_dataset       
            
    @property
    def point_utils(self):
        '''
        Getter method for point_utils property. Will set self._point_utils if None
        '''
        if not self._point_utils:
            self._point_utils = NetCDFPointUtils(self.netcdf_dataset, 
                                                 enable_disk_cache=False, 
                                                 enable_memory_cache=True,
                                                 debug=False #settings['global_settings']['debug']
                                                 )
        return self._point_utils
    
    @property
    def line_utils(self):
        '''
        Getter method for line_utils property. Will set self._line_utils if None
        '''
        if not self._line_utils:
            self._line_utils = NetCDFLineUtils(self.netcdf_dataset, 
                                               enable_disk_cache=False,
                                               enable_memory_cache=True,
                                               debug=False #settings['global_settings']['debug']
                                               )
        self._point_utils = self._line_utils # NetCDFLineUtils is a subclass of NetCDFPointUtils
        return self._line_utils
        
    @property
    def grid_utils(self):
        '''
        Getter method for grid_utils property. Will set self._grid_utils if None
        '''
        if not self._grid_utils:
            self._grid_utils = NetCDFGridUtils(self.netcdf_dataset,
                                               debug=False #settings['global_settings']['debug']
                                               )
        return self._grid_utils
    
    @property
    def netcdf_path(self):
        '''
        Getter method for netcdf_path property.
        '''
        return self._netcdf_path
    
    @netcdf_path.setter
    def netcdf_path(self, netcdf_path):
        '''
        Setter method for netcdf_path property. Will reset any dependent properties.
        '''
        if self._netcdf_path:
            if self._netcdf_dataset:
                logger.debug('Closing {}'.format(self._netcdf_path))
                self._netcdf_dataset.close()
                self._netcdf_dataset = None
            if self._point_utils:
                # del self._point_utils
                logger.debug('Setting self._point_utils = None')
                self._point_utils = None
            if self._line_utils:
                # del self._line_utils
                logger.debug('Setting self._line_utils = None')
                self._line_utils = None
            if self._grid_utils:
                # del self._grid_utils
                logger.debug('Setting self._grid_utils = None')
                self._grid_utils = None
                
        self._netcdf_path = str(netcdf_path).strip()
                
        self.netcdf_basename = os.path.basename(netcdf_path)
                
    @property
    def kml_string(self):
        '''
        Getter function to return KML string from self.dataset_type_folder or empty folder
        '''
        return str(self.dataset_type_folder) 

    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, debug_value):
        if self._debug != debug_value or self._debug is None:
            self._debug = debug_value
            
            if self._debug:
                logger.setLevel(logging.DEBUG)
                logging.getLogger(self.__module__).setLevel(logging.DEBUG)
            else:
                logger.setLevel(logging.INFO)
                logging.getLogger(self.__module__).setLevel(logging.INFO)
                
        logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))
        logging.getLogger(self.__module__).debug('Logger {} set to level {}'.format(self.__module__, logger.level))
       
        
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
#         dataset_kml = kml.newfolder(name="Ground Gravity Survey Observations")
#         netcdf2kml_obj = NetCDF2kmlConverter(survey)
#         netcdf2kml_obj.netcdf_dataset = netCDF4.Dataset(netcdf2kml_obj.netcdf_path)
#         netcdf2kml_obj.npu = NetCDFPointUtils(netcdf2kml_obj.netcdf_dataset)
# 
#         survey_region = netcdf2kml_obj.build_region()
#         polygon_region = netcdf2kml_obj.build_region(0, 100)
# 
#         dataset_kml = kml.newfolder(name="polygon")
#         dataset_kml = netcdf2kml_obj.build_polygon(dataset_kml)
#         dataset_kml.region = polygon_region  # insert built point region into point folder
# 
#         dataset_kml.region = survey_region  # insert built point region into point folder
# 
#         if netcdf2kml_obj.npu.point_count > 2:
# 
#             survey_points_folder = netcdf2kml_obj.build_points(dataset_kml,
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
