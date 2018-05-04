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
from geophys_utils.netcdf_converter import NetCDFConverter, NetCDFVariable
import sys
import re
from datetime import datetime
import yaml
import os
import logging

# # Create the Logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Create the console handler and set logging level
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
# Create a formatter for log messages
logger_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# Add the Formatter to the Handler
console_handler.setFormatter(logger_formatter)
# Add the Handler to the Logger
logger.addHandler(console_handler)


class Grav2NetCDFConverter(NetCDFConverter):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''

    gravity_metadata_list = [
        # 'ENO',
        'SURVEYID',
        'SURVEYNAME',
        'COUNTRYID',
        'STATEGROUP',
        'STATIONS',
        # 'GRAVACC', - variable
        ['GRAVDATUM', 'GRAVDATUMS'],  # TODO always 'B'? Australian Absulte Gravity Datum 2007 (AAGD07)
        # 'GNDELEVACC', - variable
        # 'GNDELEVMETH', - variable
        # 'GNDELEVDATUM', - variable - 6 outlyers
        # 'RELIAB', variable - 5 outlyers
        'LAYOUT',  # fuller descriptions of this are somewhere.
        # 'ACCESS_CODE', filtered
        # 'ENTRYDATE',
        # 'ENTEREDBY',
        # 'LASTUPDATE',
        # 'UPDATEDBY',
        # 'GRAVACCUNITS', #always um. put in grav acc var attribute - may be null sometimes
        # 'GRAVACCMETHOD', variable
        'GNDELEVACCUNITS',  # always m maybe some as null
        # 'GNDELEVACCMETHOD', as variable
        'ELLIPSOIDHGTDATUM',  # always - always GRS80
        'ELLIPSOIDHGTMETH',
        # 'ELLIPSOIDHGTACC', # as variable
        # 'ELLIPSOIDHGTACCMETHOD',# as variable
        'ELLIPSOIDHGTACCUOM',
        'SURVEYTYPE',
        # 'DATATYPES',
        # 'UNO',
        'OPERATOR',
        'CONTRACTOR',
        'PROCESSOR',
        'CLIENT',  # nulls
        'OWNER',  # nulls
        'LEGISLATION',  # nulls
        # 'STATE',
        'PROJ_LEADER',  # nulls
        'ON_OFF',
        'STARTDATE',
        'ENDDATE',
        'VESSEL_TYPE',  # nulls
        'VESSEL',  # nulls
        'SPACEMIN',  # can add uom which is metres
        'SPACEMAX',
        # 'LOCMETHOD', - variable
        'ACCURACY',  # ???
        # 'GEODETIC_DATUM',
        'PROJECTION',  # nulls
        # 'QA_CODE',
        'RELEASEDATE',  # not needed but open for discussion
        'COMMENTS',  # not needed but open for discussion
        # 'DATA_ACTIVITY_CODE',
        # 'NLAT', already in global attributes
        # 'SLAT', already in global attributes
        # 'ELONG', already in global attributes
        # 'WLONG', already in global attributes
        # 'ANO',
        # 'QABY',
        # 'QADATE',
        # 'CONFID_UNTIL',
    ]

    try:
        logger.debug(os.path.splitext(__file__)[0] + '_settings.yml')
        settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_settings.yml'))
        logger.debug('settings' + str(settings))
    except:
        logger.debug("Yaml load fail")
        settings = {}

    def get_keys_and_values(self, table_name: str):
        """
        Retrieves all data from a specified table, converts into a dictionary, and returns as a string. Used for tables
        with the key and value information such as accuray or methodology.
        e.g. 'SUR': 'Positions determined by optical surveying methods or measured on surveyed points.'
        """
        sql_statement = 'select * from gravity.{}'.format(table_name)
        query_result = self.cursor.execute(sql_statement)
        accuracy_method_keys_and_values_dict = {}
        for s in query_result:
            # for every instance in the table, add the 1st and 2nd column as key, value in a python dict
            accuracy_method_keys_and_values_dict[s[0]] = s[1]

        # returns as string. Python dict not accepted.
        return accuracy_method_keys_and_values_dict

    def get_value_for_key(self, value_column: str, table_name: str, key_column: str,  key: str):
        """
        Retrieves all data from a specified table, converts into a dictionary, and returns as a string. Used for tables
        with the key and value information such as accuray or methodology.
        e.g. 'SUR': 'Positions determined by optical surveying methods or measured on surveyed points.'
        """
        sql_statement = "select {0} from gravity.{1} where {2} = '{3}'".format(value_column, table_name, key_column, key)
        query_result = self.cursor.execute(sql_statement)
        cleaned_target_value = str(next(query_result))
        list_of_characters_to_remove = ["\(", "\)", "\'", "\,"]
        for character in list_of_characters_to_remove:
            cleaned_target_value = re.sub(character, '', cleaned_target_value)
        return cleaned_target_value

    def __init__(self, nc_out_path, survey_id, con, sql_strings_dict_from_yaml, netcdf_format='NETCDF4'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''

        NetCDFConverter.__init__(self, nc_out_path, netcdf_format)

        self.cursor = con.cursor()
        self.survey_id = survey_id
        self.sql_strings_dict_from_yaml = sql_strings_dict_from_yaml

        self.survey_metadata = self.get_survey_metadata()

    def get_survey_metadata(self):
        """
        Retrieve all data from the gravsurveys and joined a.surveys tables for the current surveyid in the loop.
        Uses same filters as other sql queries.

        :return:
        """
        # TODO are the filters needed in the sql? It will pass this survey id if no observation data is used later on?

        print(self.sql_strings_dict_from_yaml['get_survey_metadata'].format(self.survey_id))
        formatted_sql = self.sql_strings_dict_from_yaml['get_survey_metadata'].format(self.survey_id)
        query_result = self.cursor.execute(formatted_sql)

        field_names = [field_desc[0] for field_desc in query_result.description]
        survey_row = next(query_result)

        return dict(zip(field_names, survey_row))


    # def get_survey_metadata_in_obs_table(self):
    #     columns_to_add = {'LOCCACCUOM' : None }
    #     for key, value in iter(columns_to_add.items()):
    #         # sql_statement = '''select {0} from gravity.OBSERVATIONS go
    #         #                      where gravity.GRAVSURVEYS.surveyid = {1}
    #         #                      and go.dlong is not null
    #         #                      and go.dlat is not null
    #         #                      and status = 'A'
    #         #                      and access_code = 'O'
    #         #                      and geodetic_datum = 'GDA94'
    #         #                      )'''.format(key, survey_id)
    #
    #         query_result = self.cursor.execute(sql_statement)
    #         value = next(query_result)
    #     return value


    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''

        metadata_dict = {'title': self.survey_metadata['SURVEYNAME'],
            'Conventions': "CF-1.6,ACDD-1.3",
            'Gravity_Accuracy' : self.survey_metadata['GRAVACC'], #example of how to add oracle fields to global attributes
            'featureType': "trajectory",
            'keywords': 'blah',
            'geospatial_east_min': np.min(self.nc_output_dataset.variables['Long']),
            'geospatial_east_max': np.max(self.nc_output_dataset.variables['Long']),
            'geospatial_east_units': "m",
            'geospatial_east_resolution': "point",
            'geospatial_north_min': np.min(self.nc_output_dataset.variables['Lat']),
            'geospatial_north_max': np.max(self.nc_output_dataset.variables['Lat']),
            'geospatial_north_units': "m",
            'geospatial_north_resolution': "point",
            'geospatial_vertical_min': np.min(self.nc_output_dataset.variables[('Gndelev')]), # TODO say if I use gndelev or meter height
            'geospatial_vertical_max': np.max(self.nc_output_dataset.variables[('Gndelev')]), # TODO this be min(elevation-DOI)?
            'geospatial_vertical_units': "m",
            'geospatial_vertical_resolution': "point",
            'geospatial_vertical_positive': "up",
           # 'history': 'Pulled from database at Geoscience Australia'
            'date_created': datetime.now().isoformat()
            }

        return metadata_dict

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''

        formatted_sql = self.sql_strings_dict_from_yaml['get_dimensions'].format(self.survey_id)
        print(formatted_sql)
        self.cursor.execute(formatted_sql)
        point_count = int(next(self.cursor)[0])

        dimensions = OrderedDict()
        dimensions['point'] = point_count  # number of points per survey

        return dimensions



    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects

        '''

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
                    if isinstance(metadata_attribute, list):
                        if key == metadata_attribute[0]:
                            # get_value_for_key(value_column: str, table_name: str, key_column: str,  key: str)
                            gravity_metadata[key] = str(self.get_value_for_key('DESCRIPTION', metadata_attribute[1], key, value))
            logger.debug("GA gravity metadata")
            logger.debug(gravity_metadata)

            return gravity_metadata

        def get_data(field_name_dict):
            """

            :param field_name_dict:
            :return:
            """

            # call the sql query and assign results into a python list
            formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(field_name_dict['database_field_name'], self.survey_id)
            self.cursor.execute(formatted_sql)

            variable_list = []
            for i in self.cursor:
                variable_list.append(i[0])  # getting the first index is required. Otherwise each point is within its own tuple.

            return variable_list


                # return as numpy array, with dtype specified in yaml file.
                #return np.array(variable_list, dtype=field_name_dict['dtype'])


        def convert_list_to_mapped_values(list_to_edit, mapping_dict):

            transformed_list = []

            for l in list_to_edit:
                logger.debug("value: " + str(l))
                for key5, value5, in mapping_dict.items():
                    if l == key5:
                        print(mapping_dict.get(key5))
                        transformed_list.append(mapping_dict.get(key5))
                    else:
                        pass
            print('transformed_list')
            return transformed_list

        def get_field_description(target_field):
            sql_statement = """
                SELECT COMMENTS 
                FROM ALL_COL_COMMENTS   
                WHERE TABLE_NAME = 'OBSERVATIONS' 
                AND COLUMN_NAME = '{}'""".format(target_field.upper())
            self.cursor.execute(sql_statement)
            comment = str(next(self.cursor)[0])
            return comment

        def handle_key_value_cases(field_value, value):
            key_values_tables_dict = self.get_keys_and_values(field_value.get(value))
            print("DICT")
            print(key_values_tables_dict)
            print(type(key_values_tables_dict))

            key_list = []
            # get the keys into a list

            for key, value2 in key_values_tables_dict.items():
                # if re.search('[A-Z]', key):
                print(key, value2)
                key_list.append(key)

            # convert the key list to a np array
            lookup_array = np.array(key_list)
            print("lookup array")
            print(lookup_array)
            value_array = get_data(field_value)

            # create the mapping dict/lookup table to convert variables with strings as keys.
            mapping_dict = {}
            for this in key_list:
                mapping_dict[this] = key_list.index(this)

            transformed_list = convert_list_to_mapped_values(value_array, mapping_dict)
            # ummm


            if attributes_dict['dtype'] == 'S4':
                pass
                # print(attributes_dict['comments'])

            return transformed_list

        # Begin yielding NetCDFVariables

        # crs variable creation for GDA94
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



        # survey level attributes
        # add loccaccuom
        # GNDELEVACCUOM
        # METERHGTUNITS
        # GNDELEVACCUOM
        # METERHGTERRUOM
        # GRAVACCUOM - one survey is all null
        # TCERRMETHOD some are all nulls
        # ELLIPSOIDHGTUNITS - always m
        # ELLIPSOIDHGTMETH
        # ELLIPSOIDHGTDATUM - always GRS80

        # ELLIPSOIDMETERHGTUNITS is always m - added as units attribute for ellipsoide meter hgt
        # ELLIPSOIDHGTACCUOM is always m

        yield NetCDFVariable(short_name='ga_gravity_metadata',
                              data=0,
                              dimensions=[],  # Scalar
                              fill_value=None,
                              attributes=generate_ga_metadata_dict(),
                              dtype='int8'  # Byte datatype
                              )

        # values to parse into NetCDFVariable attributes list. Once passed they become a netcdf variable attribute.
        # key_value_table is later converted to comments.
        list_of_possible_value = ['long_name', 'units', 'dtype', 'key_value_table']

        for field_name, field_value in Grav2NetCDFConverter.settings['field_names'].items():
            logger.debug('-----------------')
            logger.debug("Field Names: " + str(field_name))
            logger.debug("Field Values: " + str(field_value))
            attributes_dict = {'description': get_field_description(field_value['database_field_name'])}

            for value in list_of_possible_value:
                logger.debug("Value in list_of_possible_value: " + str(value))

                # if the field value is in the list of accepted values then add to attributes dict
                if field_value.get(value):
                    # handle the key_value_table madness
                    if value == 'key_value_table':
                        # get the transformed data if the data is in string form
                        transformed_list = handle_key_value_cases(field_value, value)
                        attributes_dict['comments'] = str(self.get_keys_and_values(field_value.get(value)))
                    # otherwise it's easy
                    else:
                        attributes_dict[value] = field_value[value]
                else:
                    logger.debug(str(field_name) + ' is not set as an accepted attribute.')
                    logger.debug('attributes_dict' + str(attributes_dict))
            logger.debug('attributes_dict' + str(attributes_dict))
            #print('attributes_dict' + str(attributes_dict))
            variable_data = np.array(get_data(field_value), dtype=field_value['dtype'])

            if variable_data is not None or transformed_list is not None:

                yield NetCDFVariable(short_name=field_value['short_name'],
                                     data=variable_data if variable_data is not None else transformed_list,
                                     dimensions=['point'],
                                     fill_value=None,
                                     attributes=attributes_dict
                                     )


def main():
    # get user input and connect to oracle
    assert len(sys.argv) >= 4, '....'
    nc_out_path = sys.argv[1]
    u_id = sys.argv[2]
    oracle_database = sys.argv[3]
    pw = sys.argv[4]
    con = cx_Oracle.connect(u_id, pw, oracle_database)
    survey_cursor = con.cursor()

    # get sql strings from yaml file
    yaml_sql_settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_sql_strings.yml'))
    sql_strings_dict = yaml_sql_settings['sql_strings_dict']
    print(sql_strings_dict)
    # execute sql to return surveys to convert to netcdf
    survey_cursor.execute(sql_strings_dict['sql_get_surveyids'])
    survey_id_list = []

    # tidy the survey id strings
    for survey_row in survey_cursor:
        tidy_sur = re.search('\d+', survey_row[0]).group()
        survey_id_list.append(tidy_sur)

    logger.debug('Survey count =', len(survey_id_list))
    logger.debug(survey_id_list)

    # Loop through he survey lists to make a netcdf file based off each one.
    for survey in survey_id_list:
        if survey == '197309':
            logger.debug("Processing for survey: " + str(survey))
            g2n = Grav2NetCDFConverter(nc_out_path + str(survey) + '.nc', survey, con, sql_strings_dict)
            g2n.convert2netcdf()
            logger.info('Finished writing netCDF file {}'.format(nc_out_path))
            logger.info('Global attributes:')
            for key, value in iter(g2n.nc_output_dataset.__dict__.items()):
                logger.info(str(key) + ": " + str(value))
            #logger.info(g2n.nc_output_dataset.__dict__)
            logger.info('Dimensions:')
            logger.info(g2n.nc_output_dataset.dimensions)
            logger.info('Variables:')
            logger.info(g2n.nc_output_dataset.variables)
            logger.info(g2n.nc_output_dataset.file_format)
            for data in g2n.nc_output_dataset.variables['Gndelevaccmeth']:
                print(data)
            #print(g2n.nc_output_dataset.variables['Nvalue'])

            del g2n
            break

if __name__ == '__main__':
    main()
