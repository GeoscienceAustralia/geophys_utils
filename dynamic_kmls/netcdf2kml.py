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
from datetime import date
import numpy as np
import os
import tempfile
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils
from dynamic_kmls import cache_image_file

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class NetCDF2kmlConverter(object):
    '''
    NetCDF2kmlConverter class definition
    '''
    def __init__(self, settings, dataset_type, request_host=None, debug=False):
        '''
        Constructor for NetCDF2kmlConverter class
        @param settings: Dataset settings as read from netcdf2kml_settings.yml settings file
        @param dataset_type: String indicating dataset type
        @param request_host: Optional host string needed for cached image URL
        @param debug: Boolean parameter used to turn debug output on/off
        '''
        # Initialise and set debug property
        self._debug = None
        self.debug = debug
        
        self.request_host = request_host
        
        self.cache_dir = os.path.join((settings['global_settings'].get('cache_root_dir') or 
                          tempfile.gettempdir()),
                          'kml_server_cache',
                          dataset_type
                          )
        os.makedirs(self.cache_dir, exist_ok=True)
        
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
                                    'grid': self.build_thumbnail_image
                                    }
        
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
        
           
    def build_region(self, 
                     dataset_metadata_dict,
                     min_lod_pixels=100, 
                     max_lod_pixels=-1, 
                     min_fade_extent=200, 
                     max_fade_extent=800
                     ):
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

        region = simplekml.Region(latlonaltbox="<north>" + str(dataset_metadata_dict['latitude_max']) + "</north>" +
                                               "<south>" + str(dataset_metadata_dict['latitude_min']) + "</south>" +
                                               "<east>" + str(dataset_metadata_dict['longitude_max']) + "</east>" +
                                               "<west>" + str(dataset_metadata_dict['longitude_min']) + "</west>",
                                  lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
                                      "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
                                      "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
                                      "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")
        return region

    def build_polygon(self, 
                      dataset_metadata_dict, 
                      bounding_box, 
                      visibility=True, 
                      parent_folder=None, 
                      polygon_name=None):
        """
        Builds a kml polygon into the parent folder. Polygon is built from convex_hull_polygon in database search result.
        @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @param parent_folder: Optional KML folder under which to build geometry for line dataset
        @param polygon_name: Optional polygon name - defaults to title
        
        @return: Dataset folder under parent folder
        """
        if parent_folder is None:
            parent_folder=self.dataset_type_folder
        if polygon_name is None:
            polygon_name = str(dataset_metadata_dict['dataset_title'])

        try:
            if dataset_metadata_dict['convex_hull_polygon']:
                polygon_bounds = [[float(ordinate)
                                   for ordinate in coord_pair.strip().split(' ')
                                   ]
                                  for coord_pair in
                                  re.search('POLYGON\(\((.*)\)\)',
                                            dataset_metadata_dict['convex_hull_polygon']
                                            ).group(1).split(',')
                                  ]
                # build the polygon based on the bounds. Also set the polygon name. It is inserted into the self.dataset_type_folder.
                polygon_kml = parent_folder.newpolygon(name=polygon_name,
                                               outerboundaryis=polygon_bounds, visibility=visibility)
                
                polygon_kml.style = self.polygon_style
    
                # Always set timestamps on polygons
                self.set_timestamps(polygon_kml, dataset_metadata_dict)

                # build the polygon description
                description_string = '<![CDATA['
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name',
                                                                                          str(dataset_metadata_dict['dataset_title']))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', str(dataset_metadata_dict['ga_survey_id']))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Start Date',
                                                                                          str(dataset_metadata_dict['start_date']))
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey End Date',
                                                                                          str(dataset_metadata_dict['end_date']))
                if dataset_metadata_dict['dataset_link']:
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
                    dataset_metadata_dict['dataset_link']))
                description_string = description_string + ']]>'
                polygon_kml.description = description_string
    
                return polygon_kml
        
        except Exception as e:
            logger.debug('Unable to display polygon "{}": {}'.format(dataset_metadata_dict['convex_hull_polygon'], e))
    

    def build_lines(self, dataset_metadata_dict, bounding_box, visibility=True):
        '''
        Builds a set of kml linestrings into the parent folder. Linestrings are built from dynamically subsampled points in line.
        @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        '''
        line_utils = NetCDFLineUtils(dataset_metadata_dict['netcdf_path'], 
                                               enable_disk_cache=self.cache_coordinates,
                                               enable_memory_cache=True,
                                               cache_dir=self.cache_dir,
                                               debug=self.debug
                                               )        
        # Compute segment length as a proportion of the height of bounding box
        subsampling_distance = (bounding_box[3] - bounding_box[1]) / self.line_segments_across_bbox
        
        if self.height_variable:
            height_variable = self.height_variable # e.g. 'lidar'
        else:
            height_variable = [] # Empty string to return no variables, just 'coordinates'
        
        #dataset_folder_kml.style = self.dataset_type_folder.style
        #dataset_folder_kml.style = self.line_style
    
        dataset_folder_kml = None
        
        visible_line_count = 0
        for line_number, line_data in line_utils.get_lines(line_numbers=None, 
                                                                variables=height_variable, 
                                                                bounds=bounding_box,
                                                                subsampling_distance=subsampling_distance
                                                                ):
            #logger.debug("line_number: {}".format(line_number))
            #logger.debug("line_data: {}".format(line_data))
            points_in_subset = len(line_data['coordinates'])
            if points_in_subset:
                visible_line_count += 1
                
                # Create dataset_folder_kml as soon as we have one line, leav as None if no lines found
                if not dataset_folder_kml:
                    dataset_folder_kml = self.dataset_type_folder.newfolder(name=dataset_metadata_dict['dataset_title'], visibility=visibility)
                    if self.timestamp_detail_view:
                        # Enable timestamps on lines
                        self.set_timestamps(dataset_folder_kml, dataset_metadata_dict)
        
                line_string = dataset_folder_kml.newlinestring(name=str("Line number: {}".format(line_number)))
                
                #line_string.style = dataset_folder_kml.style

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
                description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ', str(dataset_metadata_dict['dataset_title']))
                if dataset_metadata_dict['dataset_link']:
                    description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
                    dataset_metadata_dict['dataset_link']))
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
           
        line_utils.netcdf_dataset.close() # Explicitly close netCDF file
        
        if visible_line_count:
            dataset_folder_kml.name = dataset_folder_kml.name + ' ({} lines in view)'.format(visible_line_count)
            
        return dataset_folder_kml
    

    def build_points(self, dataset_metadata_dict, bounding_box, visibility=True):
        """
        Builds all points for a survey. Including building the containing folder, setting time stamps, setting the
         style, and setting the description html to pop up when the point is selected.
        @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        """        
        point_utils = NetCDFPointUtils(dataset_metadata_dict['netcdf_path'], 
                                       enable_disk_cache=self.cache_coordinates, 
                                       enable_memory_cache=True,
                                       cache_dir=self.cache_dir,
                                       debug=self.debug
                                       )
        
        spatial_mask = point_utils.get_spatial_mask(bounding_box)
        #logger.debug('spatial_mask: {}'.format(spatial_mask))
        if not np.any(spatial_mask):
            logger.debug('No points in view')
            return
        
        dataset_folder_kml = self.dataset_type_folder.newfolder(name=dataset_metadata_dict['dataset_title'], visibility=visibility)
        
        dataset_folder_kml.style = self.point_style

        if self.timestamp_detail_view:
            # Enable timestamps on points
            self.set_timestamps(dataset_folder_kml, dataset_metadata_dict)
            
        point_data_generator = point_utils.all_point_data_generator(self.point_field_list, spatial_mask)
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
            point_kml = dataset_folder_kml.newpoint(name="Observation no. " + str(point_data['obsno']),
                                                    coords=[(point_data['longitude'], point_data['latitude'])],
                                                    visibility=visibility)

            point_kml.style = dataset_folder_kml.style
            
            description_string = self.build_html_description_string(dataset_metadata_dict, variable_attributes, point_data)
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
                
        point_utils.netcdf_dataset.close() # Explicitly close netCDF file

        dataset_folder_kml.region = self.build_region(dataset_metadata_dict, 100, -1, 200, 800)
        
        if visible_point_count:
            dataset_folder_kml.name = dataset_folder_kml.name + ' ({} points in view)'.format(visible_point_count)
                        
        return dataset_folder_kml


    def build_thumbnail_image(self, dataset_metadata_dict, bounding_box, visibility=True):
        """
        Builds a kml thumbnail image into the parent folder. 
        Thumbnail URL is built from OPeNDAP endpoint at this stage, but this needs to change.
        @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bounding_box: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        @return: Dataset folder under parent folder
        """
        logger.debug("Building WMS thumbnail...")

        #=======================================================================
        # grid_utils = NetCDFGridUtils(dataset_metadata_dict['netcdf_path'],
        #                              debug=self.debug
        #                              )        
        #=======================================================================
        
        dataset_folder_kml = self.dataset_type_folder.newfolder(name=dataset_metadata_dict['dataset_title'], visibility=True)

        transparent_polygon = self.build_polygon(dataset_metadata_dict,
                                                 bounding_box, visibility=True, 
                                                 parent_folder=dataset_folder_kml,
                                                 polygon_name=dataset_folder_kml.name
                                                 )
        logger.debug('transparent_polygon: {}'.format(transparent_polygon))
        #transparent_polygon.color =
        transparent_polygon.style.polystyle.color = '03000000' # 99% transparent black
        transparent_polygon.style.polystyle.outline = 0  # remove the outline
        #transparent_polygon.style.linestyle.color = '80f8f8ff' # 50% transparent white

        try:
            logger.debug("Dataset WEST extent: {}".format(dataset_metadata_dict['longitude_min']))
            logger.debug("BBOX WEST extent: {}".format(bounding_box[0]))
            logger.debug("Dataset EAST extent: {}".format(dataset_metadata_dict['longitude_max']))
            logger.debug("BBOX EAST extent: {}".format(bounding_box[2]))
            logger.debug("Dataset SOUTH extent: {}".format(dataset_metadata_dict['latitude_min']))
            logger.debug("BBOX SOUTH extent: {}".format(bounding_box[1]))
            logger.debug("Dataset NORTH extent: {}".format(dataset_metadata_dict['latitude_max']))
            logger.debug("BBOX NORTH extent: {}".format(bounding_box[3]))

            wms_url = dataset_metadata_dict['distribution_url'].replace('/dodsC/', '/wms/') #TODO: Replace this hack

            if self.cache_images and self.request_host:
                # Retrieve image for entire dataset
                north = dataset_metadata_dict['latitude_max']
                south = dataset_metadata_dict['latitude_min']
                east = dataset_metadata_dict['longitude_max']
                west = dataset_metadata_dict['longitude_min']
            else:  
                # Retrieve image for portion of dataset in view bounding box            
                west = max(bounding_box[0], dataset_metadata_dict['longitude_min'])
                east = min(bounding_box[2], dataset_metadata_dict['longitude_max'])
                south = max(bounding_box[1], dataset_metadata_dict['latitude_min'])
                north = min(bounding_box[3], dataset_metadata_dict['latitude_max'])

            wms_url = wms_url + "?SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap&BBOX={0},{1},{2},{3}&CRS=EPSG:4326&WIDTH={4}&HEIGHT={5}&LAYERS={6}&STYLES=&FORMAT=image/png" \
                      "&DPI=120&MAP_RESOLUTION=120&FORMAT_OPTIONS=dpi:120&TRANSPARENT=TRUE" \
                      "&COLORSCALERANGE={7}%2C{8}&NUMCOLORBANDS=127".format(south, 
                                                                             west, 
                                                                             north, 
                                                                             east, 
                                                                             int((east - west) / self.wms_pixel_size), 
                                                                             int((north - south) / self.wms_pixel_size), 
                                                                             self.wms_layer_name,
                                                                             self.wms_color_range[0],
                                                                             self.wms_color_range[1]
                                                                             )
            logger.debug('wms_url: {}'.format(wms_url))

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
            
            if self.cache_images and self.request_host:
                # Cache image and mModify URL for cached image file
                wms_url = 'http://{}{}'.format(self.request_host,
                    cache_image_file(dataset_type=self.dataset_type, 
                                     image_basename=os.path.splitext(dataset_metadata_dict['netcdf_basename'])[0]+'.png', 
                                     image_source_url=wms_url)
                    )
                logger.debug('wms_url: {}'.format(wms_url))
            logger.debug('wms_url: {}'.format(wms_url))

            ground_overlay_kml = dataset_folder_kml.newgroundoverlay(name="Survey Thumbnail Image")
            ground_overlay_kml.icon.href = wms_url
            
            ground_overlay_kml.latlonbox.north = dataset_metadata_dict['latitude_max']
            ground_overlay_kml.latlonbox.south = dataset_metadata_dict['latitude_min']
            ground_overlay_kml.latlonbox.east = dataset_metadata_dict['longitude_max']
            ground_overlay_kml.latlonbox.west = dataset_metadata_dict['longitude_min']
            ground_overlay_kml.color = 'aaffffff'

            logger.debug('ground_overlay_kml.latlonbox: {}'.format(ground_overlay_kml.latlonbox))
            logger.debug('ground_overlay_kml: {}'.format(ground_overlay_kml))

            if self.timestamp_detail_view:
                self.set_timestamps(ground_overlay_kml, dataset_metadata_dict)

            logger.debug('ground_overlay_kml: {}'.format(ground_overlay_kml))
            return dataset_folder_kml
        
        except Exception as e:
            logger.debug('Unable to display thumbnail "{}": {}'.format(wms_url, e))
            pass

    

    def build_html_description_string(self, dataset_metadata_dict, variable_attributes, point_data):
        """
            Helper function to build the description string automatically based on how many fields are in the field list.
            It will take into account if a unit of measure needs to be specified or not.
            @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
            @param variable_attributes: a list of lists containing variable attributes. Matches self.point_field_list.
            @param point_data: A dict of the actual data of each point
            
            @return: return the description string html ready to be attached to a kml point.
        """

        logger.debug(variable_attributes)

        description_string = '<![CDATA['
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey Name', dataset_metadata_dict['dataset_title'])
        description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Survey ID', dataset_metadata_dict['ga_survey_id'])

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

        if dataset_metadata_dict['dataset_link']:
            description_string = description_string + '<p><b>{0}: </b>{1}</p>'.format('Link to dataset', str(
            dataset_metadata_dict['dataset_link']))
            
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
    #     logger.debug("Point Count: {}".format(point_utils.point_count))
    #     logger.debug("Area: {}".format((self.longitude_max - self.longitude_min) * (self.latitude_max - self.latitude_min)))
    #     logger.debug("Density: {}".format(((self.longitude_max - self.longitude_min) * (
    #     self.latitude_max - self.latitude_min)) / point_utils.point_count * 1000))
    #===========================================================================


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


    def set_timestamps(self, kml_entity, dataset_metadata_dict):
        '''
        Function to set timestamps on specified kml_entity
        @param kml_entity: KML entity on which to set timestamps
        @param dataset_metadata_dict_list: List of dicts containing dataset metadata as returned by DatasetMetadataCache.search_dataset_distributions function
        '''
        # Add timestamp
        # if dataset_metadata_dict['start_date'] is not None and dataset_metadata_dict['end_date'] is not None:
        try:
            assert dataset_metadata_dict['start_date'] > date(1900, 1, 1), 'Start date {} is less than 1900-01-01'.format(
                dataset_metadata_dict['start_date'])
            assert dataset_metadata_dict['end_date'] < date(2020, 1, 1), 'End date {} is greater than 2020-01-01'.format(dataset_metadata_dict['end_date'])
            assert dataset_metadata_dict['end_date'] > date(1900, 1, 1), 'End date {} is less than 1900-01-01'.format(dataset_metadata_dict['end_date'])
            assert dataset_metadata_dict['start_date'] < date(2020, 1, 1), 'Start date {} is greater than 2020-01-01'.format(
                dataset_metadata_dict['start_date'])
            kml_entity.timespan.begin = str(dataset_metadata_dict['start_date'])
            kml_entity.timespan.end = str(dataset_metadata_dict['end_date'])

        except:  # if survey does not contain start/end date information, use survey id as start/end date year.
            try:
                if dataset_metadata_dict['ga_survey_id']:
                    survey_year = int(re.match('^[0-9]{4}', str(dataset_metadata_dict['ga_survey_id'])).group())
                    assert survey_year > 1900 and survey_year < 2020, 'survey_year <= 1900 or survey_year >= 2020'
    
                    kml_entity.timespan.begin = str(survey_year) + "-06-01"
                    kml_entity.timespan.end = str(survey_year) + "-07-01"
            except Exception as e:
                logger.debug('Unable to parse year from survey ID: {}'.format(e))
                try:
                    # Look for year string and optional follow-up year string
                    year_match = re.match('.*((19|20)\d{2})((/|-)(\d{1,4}))?(\D*)', dataset_metadata_dict['dataset_title'])
                    if year_match:
                        survey_year1 = int(year_match.group(1))
                        assert survey_year1 > 1900 and survey_year1 < 2020, '{} has survey_year1 ({}) <= 1900 or survey_year1 >= 2020'.format(dataset_metadata_dict['dataset_title'], survey_year1)
                        
                        if year_match.group(5): # Survey split over two or more years
                            try:
                                survey_year2 = int(year_match.group(1)[0:5-len(year_match.group(5))]+year_match.group(5)) # Substitute year2 characters at end of year one string
                                assert survey_year2 > 1900 and survey_year2 < 2020, '{} has survey_year2 ({}) <= 1900 or survey_year2 >= 2020'.format(dataset_metadata_dict['dataset_title'], survey_year2)
                                
                                kml_entity.timespan.begin = str(survey_year1) + "-11-01"
                                kml_entity.timespan.end = str(survey_year2) + "-02-01"
                                return
                            except Exception as e:
                                logger.debug('Unable to parse year 2 from end of dataset_title: {}'.format(e))
                                
                        kml_entity.timespan.begin = str(survey_year1) + "-06-01"
                        kml_entity.timespan.end = str(survey_year1) + "-08-01"
                except Exception as e:
                    logger.debug('Unable to parse year 1 from end of dataset_title: {}'.format(e))


    def build_dataset_kml(self, kml_format, dataset_metadata_dict, bbox_list, visibility=True):
        '''
        Function to build and return KML of specified format for specified dataset
        @param kml_format: Format of KML required. Must be in ['polygon', 'point', 'line', 'grid']
        @param dataset_metadata_dict: Dict containing dataset metadata largely as returned by DatasetMetadataCache.search_dataset_distributions function
        @param bbox_list: Bounding box specified as [<xmin>, <ymin>, <xmax>, <ymax>] list
        @param visibilty: Boolean flag indicating whether dataset geometry should be visible
        '''
        build_kml_function = self.build_kml_functions.get(kml_format)
        assert build_kml_function, 'Invalid dataset form "{}". Must be in {}'.format(self.dataset_format, 
                                                                                     list(self.build_kml_functions.keys()))
        
        logger.debug('Processing {}s for dataset {}'.format(self.dataset_format, dataset_metadata_dict['netcdf_path']))
        return build_kml_function(dataset_metadata_dict, bbox_list, visibility)

        
        
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

        dataset_count = 0
        for dataset_metadata_dict in dataset_metadata_dict_list: 
            # N.B: Could determine visibility from data here
            if self.build_dataset_kml(kml_format, dataset_metadata_dict, bbox_list, visibility):   
                dataset_count += 1
        
        if dataset_count:
            self.dataset_type_folder.name = '{} {} in view'.format(dataset_count, self.dataset_type_name)
            
               
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
