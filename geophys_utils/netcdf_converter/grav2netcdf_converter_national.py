#!/usr/bin/env python

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
CSV2NetCDFConverter concrete class for converting data to netCDF

Created on 28Mar.2018

@author: Andrew Turner
'''
#TODO update creationg date

from collections import OrderedDict
import numpy as np
import cx_Oracle
from geophys_utils.netcdf_converter import ToNetCDFConverterNational, NetCDFVariableNational
from geophys_utils import points2convex_hull
import sys
import re
from datetime import datetime
import yaml
import os
import logging
import netCDF4


# # Create the Logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Create the console handler and set logging level
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.DEBUG)
# Create a formatter for log messages
logger_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# Add the Formatter to the Handler
console_handler.setFormatter(logger_formatter)
# Add the Handler to the Logger
logger.addHandler(console_handler)


class Grav2NetCDFConverter(ToNetCDFConverterNational):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''

    gravity_metadata_list = [
        # 'ENO', not needed
        #'SURVEYID', #already in global attributes
        #'SURVEYNAME', #already in global attributes as title
        'COUNTRYID',
        'STATEGROUP',
        #'STATIONS', #number of stations? Same as netcdf point dimension.
        # 'GRAVACC', - variable
        #'GRAVDATUM' as variable attribute
        # 'GNDELEVACC', - variable
        # 'GNDELEVMETH', - variable
        # 'GNDELEVDATUM', - variable - 6 outliers
        # 'RELIAB', variable - 5 outliers
        'LAYOUT',
        # 'ACCESS_CODE', filtered
        # 'ENTRYDATE', not needed
        # 'ENTEREDBY', not needed
        # 'LASTUPDATE', not needed
        # 'UPDATEDBY', not needed
        # 'GRAVACCUNITS', #always um. In grav acc var attribute - may be null sometimes
        # 'GRAVACCMETHOD', variable
        # 'GNDELEVACCUNITS',  # always m maybe some as null. In gravlevacc var attribut
        # 'GNDELEVACCMETHOD', as variable
        # 'ELLIPSOIDHGTDATUM',  # always - always GRS80 now as variable attriubte of ellipsoidhgt
        # 'ELLIPSOIDHGTMETH', methods deemed not required for analysis
        # 'ELLIPSOIDHGTACC', # as variable
        # 'ELLIPSOIDHGTACCMETHOD',# methods deemed not required for analysis
        # 'ELLIPSOIDHGTACCUOM', # as variable attribute
        'SURVEYTYPE',
        # 'DATATYPES', not needed
        # 'UNO', not needed
        'OPERATOR',
        'CONTRACTOR',
        'PROCESSOR',
        'CLIENT',  # nulls
        'OWNER',  # nulls
        'LEGISLATION',  # nulls
        # 'STATE',
        #'PROJ_LEADER',  # nulls
        'ON_OFF', #?
        #'STARTDATE', moved to global attributes for enhanced searching
        #'ENDDATE', moved to global attributes for enhanced searching
        'VESSEL_TYPE',  # nulls
        'VESSEL',  # nulls
        'SPACEMIN',  # name is changed to SPACEMIN_METRES in the post process script
        'SPACEMAX', # name is changed to SPACEMAX_METRES in the post process script
        # 'LOCMETHOD', -not needed
        #'ACCURACY',  as point variable
        #'GEODETIC_DATUM',
        #'PROJECTION', # the data is given in the netcdf as gda94 unprojected. The values in the projection are Ellispoids
        # 'QA_CODE', not needed
        #'RELEASEDATE',  # not needed but open for discussion
        'COMMENTS',  # not needed but open for discussion
        # 'DATA_ACTIVITY_CODE',
        # 'NLAT', already in global attributes
        # 'SLAT', already in global attributes
        # 'ELONG', already in global attributes
        # 'WLONG', already in global attributes
        # 'ANO', not needed
        # 'QABY', not needed
        # 'QADATE', not needed
        # 'CONFID_UNTIL', not needed
    ]

    try:
        logger.debug(os.path.splitext(__file__)[0] + '_settings.yml')
        settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_settings.yml'))
        logger.debug('Settings' + str(settings))
    except:
        logger.debug("Yaml load fail")
        settings = {}

    def get_all_ellipsoiddatum_values_table_from_point_data(self):
        """
        Search each survey and retrieve the dinstict results for ellidpoisddatum. Use these to make a lookup table,
        following the existing pattern. Currently all values will be GRS80. However, better not to hardcode it,
        just in case.
        :return: dict "{"GRS80": "GRS80"}"
        """

    def get_keys_and_values_table(self, table_name: str):
        """
        Retrieves all data from a specified table, converts into a dictionary, and returns as a string. Used for tables
        with the key and value information such as accuracy or methodology.
        e.g. 'SUR': 'Positions determined by optical surveying methods or measured on surveyed points.'
        """
        # TODO fix this
        if table_name == "ELLIPSOIDHGTDATUM":
            #return {"GRS80": "GRS80"}
            sql_statement = self.sql_strings_dict_from_yaml['get_ellipsoidhgt_datums_lookup']
            query_result = self.cursor.execute(sql_statement)
            keys_and_values_dict = OrderedDict()
            for s in query_result:
                # for every instance in the table, add the 1st and 2nd column as key, value in a python dict
                if s[0] != None:
                    keys_and_values_dict[s[0]] = s[0]
            if len(keys_and_values_dict) > 1:
                logger.debug("More than one value for ellipsoiddatum found.")
            print("ELLIPSOIDHGTDATUM")
            print("keys_and_values_dict: {}".format(keys_and_values_dict))
            return keys_and_values_dict


        sql_statement = 'select * from gravity.{}'.format(table_name)
        query_result = self.cursor.execute(sql_statement)
        keys_and_values_dict = OrderedDict()
        for s in query_result:
            # for every instance in the table, add the 1st and 2nd column as key, value in a python dict
            keys_and_values_dict[s[0]] = s[1]

        # returns as string. Python dict not accepted.
        return keys_and_values_dict

    def get_value_for_key(self, value_column: str, table_name: str, key_column: str,  key: str):
        """
        Retrieves all data from a specified table, converts into a dictionary, and returns as a string. Used for tables
        with the key and value information such as accuracy or methodology.
        e.g. 'SUR': 'Positions determined by optical surveying methods or measured on surveyed points.'
        """

        cleaned_key = str(key)
        list_of_characters_to_remove = ["\(", "\)", "\'", "\,"]
        for character in list_of_characters_to_remove:
            cleaned_key = re.sub(character, '', cleaned_key)

        sql_statement = "select {0} from gravity.{1} where {2} = '{3}'".format(value_column, table_name, key_column, cleaned_key)
        query_result = self.cursor.execute(sql_statement)
        key_to_return = str(next(query_result))

        for character in list_of_characters_to_remove:
            key_to_return = re.sub(character, '', key_to_return)

        return key_to_return

    def __init__(self, nc_out_path, survey_id, survey_index, con, sql_strings_dict_from_yaml,
                 netcdf_format='NETCDF4', national=True):
        """
        Concrete constructor for subclass Grav2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        """

        ToNetCDFConverterNational.__init__(self, nc_out_path, national, netcdf_format)

        self.survey_id = survey_id
        self.survey_index = survey_index
        self.cursor = con.cursor()
        self.sql_strings_dict_from_yaml = sql_strings_dict_from_yaml

        self.survey_metadata = self.get_survey_metadata()
        self.elipsoid_height_datums = []
        # self.get_ellipsoid_height_datum_keys()

        formatted_sql = self.sql_strings_dict_from_yaml['get_dimensions'].format(self.survey_id)
        self.cursor.execute(formatted_sql)
        self.point_count = int(next(self.cursor)[0])

    def get_survey_metadata(self):
        """
        Retrieve all data from the gravsurveys and joined a.surveys tables for the current surveyid in the loop.
        Uses same filters as other sql queries.

        :return:
        """
        # TODO are the filters needed in the sql? It will pass this survey id if no observation data is used later on?

        formatted_sql = self.sql_strings_dict_from_yaml['get_survey_metadata'].format(self.survey_id)
        query_result = self.cursor.execute(formatted_sql)
        field_names = [field_desc[0] for field_desc in query_result.description]
        survey_row = next(query_result)

        return dict(zip(field_names, survey_row))

    def get_national_survey_metadata(self):
        """
        Retrieve all survey metadata for gravity points from the db and return as a list of strings.

        :return:
        """
        sql = self.sql_strings_dict_from_yaml['get_national_survey_metadata']
        query_result = self.cursor.execute(sql)
        field_names = [field_desc[0] for field_desc in query_result.description]
        zipped_data_dicts = [dict(zip(field_names, row)) for row in query_result.fetchall()]

        # Calculate and add the time duration to each dict and then convert the start and end dates to strings
        for survey_metadata_dict in zipped_data_dicts:
            survey_metadata_dict['DURATION'] = str(survey_metadata_dict['ENDDATE'] - survey_metadata_dict['STARTDATE'])
            survey_metadata_dict['STARTDATE'] = str(survey_metadata_dict['STARTDATE'])
            survey_metadata_dict['ENDDATE'] = str(survey_metadata_dict['ENDDATE'])

        # Convert the data dicts to strings as netCDF variables do not support dicts
        zipped_data_strings = [str(data_dict) for data_dict in zipped_data_dicts]

        return zipped_data_strings

    def append_data2netcdf(self):
        logger.debug("Running append_data2netcdf function on survey: " + str(self.survey_id))
        # Get the survey metadata index variable
        smi_var = self.nc_output_dataset.variables['survey_metadata_index']
        smi_var[:] = np.append(smi_var[:], [self.survey_index] * self.point_count)
        #logger.debug(smi_var)

        for field_name, field_value in Grav2NetCDFConverter.settings['field_names'].items():
            if field_name in ['Obsno', 'Stationno', 'Lat', 'Long', 'Locacc', 'Freeair', 'Bouguer', 'Grav', 'Gravacc',
                              'Gndelev', 'Gndelevacc', 'Insthgt', 'Insthgterr', 'Ellipsoidinsthgt',
                              'Ellipsoidinsthgterr', 'Ellipsoidhgt', 'Ellipsoidhgtacc', 'Tc', 'Tcdensity', 'Tcerr']:
                if field_name in ['Freeair', 'Bouguer']:
                    formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(
                        field_value['database_field_name'],
                        field_value['fill_value'], self.survey_id)
                else:
                    formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(
                        'o1.' + field_value['database_field_name'],
                        field_value['fill_value'],
                        self.survey_id)
                try:
                    self.cursor.execute(formatted_sql)
                except:
                    logger.debug(formatted_sql)
                    raise

                print("cursor")
                print(type(self.cursor))
                data_list = [x[0] for x in self.cursor]  # get the first index. Otherwise each point is within its own tuple.

                var = self.nc_output_dataset.variables[(field_value.get('standard_name') or field_value['short_name']).lower()]
                var[:] = np.append(var[:-self.point_count], data_list)

            if field_name in ['Stattype', 'Locmethod', 'Locaccmethod', 'Gravmeth', 'Gravdatum', 'Gravaccmeth',
                              'Gndelevtype', 'Gndelevdatum', 'Gndelevmeth', 'Gndelevaccmethod', 'Insthgtmeth',
                              'Insthgterrmeth', 'Ellipsoidinsthgterrmethod', 'Ellipsoidhgtmeth',
                              'Ellipsoidhgtaccmethod', 'Ellipsoiddatum', 'Tcmeth', 'Gridflag', 'Reliab']:

                key_value_dict = self.get_keys_and_values_table(field_value.get('lookup_table'))
                logger.debug(key_value_dict)

                lookup_key_list = [lookup_key for lookup_key in key_value_dict.keys()]
                logger.debug(lookup_key_list)

                lookup_dict = {lookup_key: lookup_key_list.index(lookup_key) for lookup_key in lookup_key_list}
                logger.debug(lookup_dict)

                formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(
                    'o1.' + field_value['database_field_name'],
                    field_value['fill_value'],
                    self.survey_id)

                try:
                    self.cursor.execute(formatted_sql)
                except:
                    logger.debug(formatted_sql)
                    raise

                print("cursor")
                print(type(self.cursor))
                data_list = [x[0] for x in self.cursor]  # get the first index. Otherwise each point is within its own tuple.
                logger.debug(data_list)

                # transform the data_list into the mapped value.
                transformed_data_list = [lookup_dict.get(lookup_key) for lookup_key in data_list]
                logger.debug(transformed_data_list)
                logger.debug(len(transformed_data_list))

                for i, value in enumerate(transformed_data_list):
                    if value is None:
                        transformed_data_list[i] = int(field_value['fill_value'])

                logger.debug(transformed_data_list)

                '''converted_attributes_dict = {lookup_table_dict[key]: value
                                             for key, value in lookup_table_dict.items()}'''

                var = self.nc_output_dataset.variables[((field_value.get('standard_name') or field_value['short_name']) + '_index').lower()]
                var[:] = np.append(var[:-self.point_count], transformed_data_list)

                '''converted_data_list, converted_key_value_dict = handle_key_value_cases(field_yml_settings_dict,
                                                                                       self.get_keys_and_values_table(
                                                                                           field_yml_settings_dict.get(
                                                                                               'lookup_table')))'''


    def get_ellipsoid_height_datum_keys(self): # only one currently, GRS80
        formatted_sql = self.sql_strings_dict_from_yaml['get_ellipsoidhgt_datums_lookup']
        query_result = self.cursor.execute(formatted_sql)
        keys_and_values_dict = OrderedDict()
        for s in query_result:
            # for every instance in the table, add the 1st and 2nd column as key, value in a python dict
            keys_and_values_dict[s[0]] = s[1]
        print(keys_and_values_dict)
        self.elipsoid_height_datums = [datum[0] for datum in survey_row if datum[0] is not None]
        #return (self.elipsoid_height_datums)

    def get_survey_wide_value_from_obs_table(self, field):
        """
        Helper function to retrieve a survey wide value from the observations table. The returning value is tested
        to be the only possible value (or null) within that survey.
        :param field: The target column in the observations table.
        :return: The first value of the specified field of the observations table.
        """
        formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format('o1.'+field, "null", self.survey_id)
        formatted_sql = formatted_sql.replace('select', 'select distinct', 1)  # Only retrieve distinct results
        formatted_sql = re.sub('order by .*$', '', formatted_sql) # Don't bother sorting
        query_result = self.cursor.execute(formatted_sql)
        value = None

        for result in query_result:
            logger.debug('value: {}, result: {}'.format(value, result))
            
            if value is None:
                value = result[0]
            
            assert value is None or result[0] == value or result[0] is None, 'Variant value found in survey-wide column {}'.format(field)
        return value

    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''
        start_date_sql = self.sql_strings_dict_from_yaml['get_national_survey_metadata_order_by_startdate']
        start_date_query_result = self.cursor.execute(start_date_sql)
        field_names = [field_desc[0] for field_desc in start_date_query_result.description]
        start_date_zipped_data = dict(zip(field_names, start_date_query_result.fetchone()))
        earliest_start_date = start_date_zipped_data['STARTDATE']

        end_date_sql = self.sql_strings_dict_from_yaml['get_national_survey_metadata_order_by_enddate_desc']
        end_date_query_result = self.cursor.execute(end_date_sql)
        field_names = [field_desc[0] for field_desc in end_date_query_result.description]
        end_date_zipped_data = dict(zip(field_names, end_date_query_result.fetchone()))
        latest_end_date = end_date_zipped_data['ENDDATE']

        logger.debug(start_date_zipped_data)
        logger.debug(earliest_start_date)
        logger.debug(end_date_zipped_data)
        logger.debug(latest_end_date)
        logger.debug(latest_end_date - earliest_start_date)

        metadata_dict = {
            'title': "Australian National Ground Gravity Compilation 2024",
            'Conventions': "CF-1.6,ACDD-1.3",
            'keywords': "points, gravity, geophysical survey, Earth sciences, geophysics, geoscientific Information",
            'geospatial_lon_min': 111,
            'geospatial_lon_max': 156,
            'geospatial_lon_units': "degrees_east",
            'geospatial_lon_resolution': "point",
            'geospatial_lat_min': -6,
            'geospatial_lat_max': 45,
            'geospatial_lat_units': "degrees_north",
            'geospatial_lat_resolution': "point",
            'history': "Pulled from point gravity database at Geoscience Australia",
            'summary': "Compilation of ground gravity surveys acquired in Australia. Data acquired from State and "
                       "National Geological Surveys, Academics, and private companies. Data has gone through Quality "
                       "Control processes. Station spacing ranges from 11km to less than 1km. The accuracy of the data "
                       "varies generally with date of acquisition, later data using high precision gravity meters and "
                       "GPS for better accuracies.",
            'location_accuracy_min': np.min(self.nc_output_dataset.variables['locacc']),
            'location_accuracy_max': np.max(self.nc_output_dataset.variables['locacc']),
            'location_accuracy_units': "m",
            'elevation_accuracy_min': np.min(self.nc_output_dataset.variables['gndelevacc']),
            'elevation_accuracy_max': np.max(self.nc_output_dataset.variables['gndelevacc']),
            'elevation_accuracy_units': "m",
            'gravity_accuracy_min': np.min(self.nc_output_dataset.variables['gravacc']),
            'gravity_accuracy_max': np.max(self.nc_output_dataset.variables['gravacc']),
            'gravity_accuracy_units': "um/s^2",
            'time_coverage_start': str(earliest_start_date),
            'time_coverage_end': str(latest_end_date),
            'time_coverage_duration': str(latest_end_date - earliest_start_date),
            'date_created': datetime.now().isoformat(),
            'institution': 'Geoscience Australia',
            'source': 'ground observation',
            'cdm_data_type': 'Point'
            }

        try:
            #Compute convex hull and add GML representation to metadata
            coordinates = np.array(list(zip(self.nc_output_dataset.variables['longitude'][:],
                                            self.nc_output_dataset.variables['latitude'][:]
                                            )
                                        )
                                   )
            if len(coordinates) >=3:
                convex_hull = points2convex_hull(coordinates)        
                metadata_dict['geospatial_bounds'] = 'POLYGON((' + ', '.join([' '.join(
                    ['%.4f' % ordinate for ordinate in coordinates]) for coordinates in convex_hull]) + '))'
            elif len(coordinates) == 2: # Two points - make bounding box
                bounding_box = [[min(coordinates[:,0]), min(coordinates[:,1])],
                                [max(coordinates[:,0]), min(coordinates[:,1])],
                                [max(coordinates[:,0]), max(coordinates[:,1])],
                                [min(coordinates[:,0]), max(coordinates[:,1])],
                                [min(coordinates[:,0]), min(coordinates[:,1])]
                                ]
                metadata_dict['geospatial_bounds'] = 'POLYGON((' + ', '.join([' '.join(
                    ['%.4f' % ordinate for ordinate in coordinates]) for coordinates in bounding_box]) + '))'
            elif len(coordinates) == 1: # Single point
                #TODO: Check whether this is allowable under ACDD
                metadata_dict['geospatial_bounds'] = 'POINT((' + ' '.join(
                    ['%.4f' % ordinate for ordinate in coordinates[0]]) + '))'
        except:
            logger.warning('Unable to write global attribute "geospatial_bounds"')
            
        return metadata_dict

    def set_global_attributes(self):
        '''
        Set global attributes in netCDF output file
        '''
        for attribute_name, attribute_value in iter(self.get_global_attributes().items()):
            setattr(self.nc_output_dataset, attribute_name, attribute_value or '')

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''

        # formatted_sql = self.sql_strings_dict_from_yaml['get_dimensions'].format(self.survey_id)
        # self.cursor.execute(formatted_sql)
        # point_count = int(next(self.cursor)[0])

        dimensions = OrderedDict()
        dimensions['point'] = 0  # number of points

        formatted_sql2 = self.sql_strings_dict_from_yaml['sql_get_count_surveyids'].format(self.survey_id)
        self.cursor.execute(formatted_sql2)
        survey_count = int(next(self.cursor)[0])

        dimensions['survey'] = survey_count  # number of surveys

        for field_value in Grav2NetCDFConverter.settings['field_names'].values():
            if field_value.get('lookup_table'):
                lookup_dict = self.get_keys_and_values_table(field_value['lookup_table'])
                print("LOOKUP DICT: {}".format(lookup_dict))

                new_dimension_name = field_value['short_name'].lower()
                dimensions[new_dimension_name] = len(lookup_dict)
                # print(dimensions[new_dimension_name])
            else:
                pass
        # print(dimensions['point'])
        return dimensions

    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariableNational objects

        '''

        def get_data(field_yml_settings_dict):
            """
            Call an sql query to retrieve a data list of the specified field. A different query is called for freeair
            and bouguer.
            :param field_yml_settings_dict:
            :return: data list
            """

            if field_name in ['Freeair', 'Bouguer']:
                formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(
                    field_yml_settings_dict['database_field_name'],
                    field_yml_settings_dict['fill_value'], self.survey_id)

            else:
                formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(
                    'o1.' + field_yml_settings_dict['database_field_name'],
                    field_yml_settings_dict['fill_value'],
                    self.survey_id)

            try:
                self.cursor.execute(formatted_sql)
            except:
                logger.debug(formatted_sql)
                raise

            print("cursor")
            print(type(self.cursor))
            data_list = [x[0] for x in
                         self.cursor]  # get the first index. Otherwise each point is within its own tuple.
            # for i in self.cursor:
            #     data_list.append(
            #         i[0])
            # print(data_list)
            return data_list

        def generate_ga_metadata_dict():
            gravity_metadata = {}
            for key, value in iter(self.survey_metadata.items()):
                for metadata_attribute in Grav2NetCDFConverter.gravity_metadata_list:
                    if value is not None:
                        if key == metadata_attribute:
                            if type(value) == datetime:
                                gravity_metadata[key] = value.isoformat()
                            else:
                                gravity_metadata[key] = value

            logger.debug("GA gravity metadata")
            logger.debug(gravity_metadata)

            return gravity_metadata


        def handle_key_value_cases(field_yml_settings_dict, lookup_table_dict):
            """

            :param field_yml_settings_dict: field settings as written in the yml file e.g. {'short_name': 'Ellipsoidhgt', 'dtype':
            'float32', 'database_field_name': 'ELLIPSOIDHGT', 'long_name': 'Ellipsoid Height', 'units': 'm',
            'fill_value': -99999.9, 'datum': 'ELLIPSOIDHGTDATUM'}
            :param lookup_table_dict: Dict of key and values pulled from oracle for tables such as accuracy and
            methodology.
            :return:
            """

            logger.debug('- - - - - - - - - - - - - - - -')
            logger.debug('handle_key_value_cases() with field value: ' + str(field_value) + ' and lookup_table_dict: ' + str(lookup_table_dict))

            # get the keys into a list
            lookup_key_list = [lookup_key for lookup_key in lookup_table_dict.keys()]

            # create the lookup table to convert variables with strings as keys.
            lookup_dict = {lookup_key: lookup_key_list.index(lookup_key)
                           for lookup_key in lookup_key_list}
            
            # get the array of numeric foreign key values
            field_data_array = get_data(field_yml_settings_dict)

            # transform the data_list into the mapped value.
            transformed_data_list = [lookup_dict.get(lookup_key) for lookup_key in field_data_array]
                                
            # loop through the table_key_dict and the lookup table. When a match is found add the new mapped key to
            # the existing value of the table_key_dict in a new dict
            converted_attributes_dict = {lookup_table_dict[key]: value
                              for key, value in lookup_table_dict.items()}
            
            #===================================================================
            # converted_dict = {}
            # for keys, values in lookup_table_dict.items():
            #     for map_key, map_value in lookup_dict.items():
            #         if keys == map_key:
            #             converted_dict[map_value] = lookup_table_dict[keys]
            #===================================================================

            return transformed_data_list, converted_attributes_dict

        def handle_key_value_for_ellipsoid_datum(field_yml_settings_dict):
            """

            :param field_yml_settings_dict: field settings as written in the yml file e.g. {'short_name': 'Ellipsoidhgt', 'dtype':
            'float32', 'database_field_name': 'ELLIPSOIDHGT', 'long_name': 'Ellipsoid Height', 'units': 'm',
            'fill_value': -99999.9, 'datum': 'ELLIPSOIDHGTDATUM'}
            :param lookup_table_dict: Dict of key and values pulled from oracle for tables such as accuracy and
            methodology.
            :return:
            """

            logger.debug('- - - - - - - - - - - - - - - -')

            lookup_table_dict = {0: "GRS80"}
            # get the keys into a list
            lookup_key_list = [lookup_key for lookup_key in lookup_table_dict.keys()]

            # create the lookup table to convert variables with strings as keys.
            lookup_dict = {lookup_key: lookup_key_list.index(lookup_key)
                           for lookup_key in lookup_key_list}

            # get the array of numeric foreign key values
            field_data_array = get_data(field_yml_settings_dict)



            for i, value in enumerate(field_data_array):
                    # print("i: {}".format(i))
                    # print("value: {}".format(value))
                if value is "GRS80":
                    field_data_array[i] = 0
                else:
                    field_data_array[i] = int(field_yml_settings_dict['fill_value'])

            #transformed_data_list
            # transform the data_list into the mapped value.
          #  transformed_data_list = [lookup_dict.get(lookup_key) for lookup_key in field_data_array]
            # transformed_data_list = [lookup_dict.get(lookup_key) if lookup_key != field_yml_settings_dict['fill_value'] else field_yml_settings_dict['fill_value'] for lookup_key in field_data_array]
           # print("transformed_data_list: {}".format(transformed_data_list))

            # loop through the table_key_dict and the lookup table. When a match is found add the new mapped key to
            # the existing value of the table_key_dict in a new dict

            converted_attributes_dict = {lookup_table_dict[key]: value
                                         for key, value in lookup_table_dict.items()}
            print("converted_attributes_dict: {}".format(converted_attributes_dict))
            # ===================================================================
            # converted_dict = {}
            # for keys, values in lookup_table_dict.items():
            #     for map_key, map_value in lookup_dict.items():
            #         if keys == map_key:
            #             converted_dict[map_value] = lookup_table_dict[keys]
            # ===================================================================

            return field_data_array, converted_attributes_dict


        def get_field_description(target_field):
            """
            Helper function to retrieve the field description from a connected oracle database
            :param target_field:
            :return field_description:
            """
            sql_statement = self.sql_strings_dict_from_yaml['get_field_description'].format(target_field.upper())
            self.cursor.execute(sql_statement)
            field_description = str(next(self.cursor)[0])

            return field_description


        def build_attribute_dict_and_data_list_of_variables(field_yml_settings_dict):
            """
            For each field, the correct attributes are added. This is based on the grav2netcd_converter settings.
            The data is converted for values that have function as a lookup table. For each field, the relevant
            attribute dictionary, and np data array is returned.
            :param field_name: field name as written in yml file e.g. Ellipsoidhgt
            :param field_yml_settings_dict: field settings as written in the yml file e.g. {'short_name': 'Ellipsoidhgt', 'dtype':
            'float32', 'database_field_name': 'ELLIPSOIDHGT', 'long_name': 'Ellipsoid Height', 'units': 'm',
            'fill_value': -99999.9, 'datum': 'ELLIPSOIDHGTDATUM'}
            :return: for each field value, return its attribute dictionary, and data array or converted data array.
            """

            # values to parse into NetCDFVariableNational attributes list. Once passed they become a netcdf variable attribute.
            # lookup_table is later converted to comments.
            list_of_possible_value = ['long_name', 'standard_name', 'units', 'dtype', 'lookup_table', 'dem', 'datum']

            logger.debug('-----------------')
            logger.debug("Field Values: " + str(field_yml_settings_dict))

            converted_data_array = []
            attributes_dict = {}

            for value in list_of_possible_value:
                logger.debug("Value in list_of_possible_value: " + str(value))
                # if the field value is in the list of accepted values then add to attributes dict
                if field_yml_settings_dict.get(value):
                    logger.debug("Processing: " + str(value))

                    # some key values are already int8 and don't need to be converted. Thus a flag is included in the
                    if value == 'lookup_table':
                        logger.debug('Converting ' + str(value) + 'string keys to int8 with 0 as 1st index')
                        converted_data_list, converted_key_value_dict = handle_key_value_cases(field_yml_settings_dict,
                                                                                               self.get_keys_and_values_table(
                                                                                                   field_yml_settings_dict.get(
                                                                                                       'lookup_table')))

                        logger.debug('Adding converted lookup table as variable attribute...')
                        # this replaces ['comments'] values set in the previous if statement.
                        # attributes_dict['comments'] = str(converted_key_value_dict)
                        # print('this format')
                        # print('converted_data_list')
                        # print(converted_data_list)
                        print("before converted data list: {}".format(converted_data_list))
                        print(converted_data_list[0])
                        # if converted_data_list[0] == None:
                        for i, value in enumerate(converted_data_list):
                            # print("i: {}".format(i))
                            # print("value: {}".format(value))
                            if value is None:
                                converted_data_list[i] = int(field_yml_settings_dict['fill_value'])
                            else:
                                converted_data_list[i] = value
                        print("converted data list: {}".format(converted_data_list))

                        converted_data_array = np.array(converted_data_list, field_yml_settings_dict['dtype'])

                    # for the one case where a column in the observation table (tcdem) needs to be added as the
                    # attribute of varaible in the netcdf file.
                    if value == 'dem' or value == 'datum':
                        # the grav datum needs to be converted from its key value
                        if field_yml_settings_dict.get('short_name') == 'Grav':
                            gravdatum_key = self.get_survey_wide_value_from_obs_table(
                                field_yml_settings_dict.get(value))
                            attributes_dict[value] = self.get_value_for_key("DESCRIPTION", "GRAVDATUMS", "GRAVDATUM",
                                                                            gravdatum_key)
                        # while TCDEM and ELLIPSOIDHGTDATUM do not
                        else:
                            attributes_dict[value] = self.get_survey_wide_value_from_obs_table(
                                field_yml_settings_dict.get(value))
                            # if None is returned then remove the attribute
                            if attributes_dict[value] is None:
                                attributes_dict.pop(value)
                            else:
                                pass

                    # for all other values, simply add them to attributes_dict
                    else:
                        # TODO fix this
                        try:
                            attributes_dict[value] = field_yml_settings_dict[value]
                            logger.debug('attributes_dict["{}"] = {}'.format(value, field_yml_settings_dict[value]))
                        except:
                            print("no go")
                # if the value isn't in the list of accepted attributes
                else:
                    logger.debug(str(value) + ' is not found in yaml config or is not set as an accepted attribute.')

            logger.debug('Attributes_dict' + str(attributes_dict))

            # if the data array contained a lookup and was converted, return it and the attribute dict.
            if len(converted_data_array) > 0:
                return converted_data_array, attributes_dict

            # else get the non converted data and return it in an numpy array and the and the attribute dict too
            else:
                data_array = np.array(get_data(field_yml_settings_dict), dtype=field_yml_settings_dict['dtype'])
                return data_array, attributes_dict

        # ------------------------------------------------------------------------------------
        # Begin yielding NetCDFVariables
        # ------------------------------------------------------------------------------------
        
        # ---------------------------------------------------------------------------
        # crs variable creation for GDA94
        # ---------------------------------------------------------------------------
        yield self.build_crs_variable('''\
        GEOGCS["GDA94",
            DATUM["Geocentric_Datum_of_Australia_1994",
                SPHEROID["GRS 1980",6378137,298.257222101,
                    AUTHORITY["EPSG","7019"]],
                TOWGS84[0,0,0,0,0,0,0],
                AUTHORITY["EPSG","6283"]],
            PRIMEM["Greenwich",0,
                AUTHORITY["EPSG","8901"]],
            UNIT["degree",0.0174532925199433,
                AUTHORITY["EPSG","9122"]],
            AUTHORITY["EPSG","4283"]]
        '''
                                          )
        # ---------------------------------------------------------------------------
        # non acc convention survey level metadata grouped into one variable
        # ---------------------------------------------------------------------------

        '''yield NetCDFVariable(short_name='ga_gravity_metadata',
                             data=0,
                             dimensions=[],  # Scalar
                             fill_value=None,
                             attributes=generate_ga_metadata_dict(),
                             dtype='int8'  # Byte datatype
                             )'''

        # ---------------------------------------------------------------------------
        # survey level metadata for all surveys
        # ---------------------------------------------------------------------------
        # TODO: Revisit attributes
        yield NetCDFVariableNational(short_name='survey_metadata',
                             data=np.array(self.get_national_survey_metadata()),
                             dimensions=['survey'],
                             fill_value=int('-99'),
                             attributes={'long_name': 'Survey Metadata'},
                             national=True
                             )

        yield NetCDFVariableNational(short_name='survey_metadata_index',
                             data=np.array([self.survey_index] * self.point_count, dtype='int8'),
                             dimensions=['point'],
                             fill_value=int('-99'),
                             attributes={'long_name': 'zero-based index of value in survey_metadata',
                                         'lookup': 'survey_metadata'
                                         },
                             national=True
                             )

        # ---------------------------------------------------------------------------
        # The point dimension variables and their assocciated lookup table variables
        # ---------------------------------------------------------------------------
        
        # Loop through the defined variables in the yaml config and construct as netcdf variables.
        for field_name, field_value in Grav2NetCDFConverter.settings['field_names'].items():
            # convert strings to int or floats for int8 and float32 to get the required data type for the fill value
            if field_value['dtype'] == 'int8':
                fill_value = int(field_value['fill_value'])
            elif field_value['dtype'] == 'float32':
                fill_value = float(field_value['fill_value'])
            else:
                fill_value = field_value['fill_value']

            data, attributes = build_attribute_dict_and_data_list_of_variables(field_value)

            if field_value.get('lookup_table'):

                # get the values from the lookup table dict and convert into a np.array
                lookup_table_dict = self.get_keys_and_values_table(field_value['lookup_table'])
                grid_value_list = [value for value in iter(lookup_table_dict.values())]
                lookup_table_array = np.array(grid_value_list)
                attributes.pop('dtype', None)
                attributes.pop('lookup_table', None)

                dim_name = field_value['short_name'].lower()

                # Yield lookup table with same name as field
                yield NetCDFVariableNational(short_name=dim_name,
                                     data=lookup_table_array,
                                     dimensions=[dim_name],
                                     fill_value=fill_value,
                                     attributes=attributes,
                                     national=True
                                     )
                
                # Yield index table with name of <field_name>_index 
                index_attributes = dict(attributes)
                index_attributes['long_name'] = "zero-based index of value in " + dim_name
                index_attributes['lookup'] = dim_name
                    
                yield NetCDFVariableNational(short_name=((field_value.get('standard_name') or field_value['short_name']) + '_index').lower(),
                                     data=data,
                                     dimensions=['point'],
                                     fill_value=fill_value,
                                     attributes=index_attributes,
                                     national=True
                                     )
                
            else: # Not a lookup field
                yield NetCDFVariableNational(short_name=(field_value.get('standard_name') or field_value['short_name']).lower(),
                                     data=data,
                                     dimensions=['point'],
                                     fill_value=fill_value,
                                     attributes=attributes,
                                     national=True
                                     )

def main():

    # get user input and connect to oracle
    assert len(sys.argv) >= 4, '....'
    nc_out_path = sys.argv[1]
    oracle_database = sys.argv[2]
    u_id = sys.argv[3]
    pw = sys.argv[4]
    con = cx_Oracle.connect(u_id, pw, oracle_database)
    survey_cursor = con.cursor()

    # Create the initial national netCDF4 file. If it already exists the data will be clobbered and the file will be
    # emptied and ready for appending data to.
    ds = netCDF4.Dataset("{}/NATIONAL_GNDGRAV.nc".format(nc_out_path), mode='w')
    ds.close()

    # get sql strings from yaml file
    yaml_sql_settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_sql_strings.yml'))
    sql_strings_dict = yaml_sql_settings['sql_strings_dict']
    # execute sql to return surveys to convert to netcdf

    survey_cursor.execute(sql_strings_dict['sql_get_surveyids'])
    
    # tidy the survey id strings
    survey_id_list = [re.search('\d+', survey_row[0]).group()
                      for survey_row in survey_cursor
                      ]

    logger.debug('Survey count = {}'.format(len(survey_id_list)))

    for survey in survey_id_list:
        # On the first iteration of the loop create the dimensions, variables and global attributes and write the data
        # for the first survey into the variables.
        survey_index = survey_id_list.index(survey)
        if survey_index == 0:

            logger.debug("Creating initial national file on survey: " + str(survey))
            g2n = Grav2NetCDFConverter("{}/NATIONAL_GNDGRAV.nc".format(nc_out_path), survey, survey_index, con, sql_strings_dict)
            g2n.convert2netcdf()
            del g2n
        elif 0 < survey_index < 4:
            g2n = Grav2NetCDFConverter("{}/NATIONAL_GNDGRAV.nc".format(nc_out_path), survey, survey_index, con, sql_strings_dict)
            g2n.append_data2netcdf()
            del g2n
        elif survey_index == 4:
            g2n = Grav2NetCDFConverter("{}/NATIONAL_GNDGRAV.nc".format(nc_out_path), survey, survey_index, con, sql_strings_dict)
            g2n.append_data2netcdf()
            g2n.set_global_attributes()
            del g2n


if __name__ == '__main__':
    main()
