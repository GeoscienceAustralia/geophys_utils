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
console_handler.setLevel(logging.DEBUG)
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

    def get_keys_and_values_table(self, table_name: str):
        """
        Retrieves all data from a specified table, converts into a dictionary, and returns as a string. Used for tables
        with the key and value information such as accuray or methodology.
        e.g. 'SUR': 'Positions determined by optical surveying methods or measured on surveyed points.'
        """
        print("TABLE_NAME: " + str(table_name))
        sql_statement = 'select * from gravity.{}'.format(table_name)
        print(sql_statement)
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
        with the key and value information such as accuracy or methodology.
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

        def handle_key_value_cases(field_value, key_values_tables_dict):
            """
            """
            logger.debug('- - - - - - - - - - - - - - - -')
            logger.debug('handle_key_value_cases() with field value: ' + str(field_value) + ' and key_value_tables_dict: ' + str(key_values_tables_dict))
            key_list = []
            # get the keys into a list
            for key, value2 in key_values_tables_dict.items():
                key_list.append(key)

            # create the mapping dict/lookup table to convert variables with strings as keys.
            mapping_dict = {}
            for this in key_list:
                mapping_dict[this] = key_list.index(this)

            # get the data
            value_array = get_data(field_value)
            # transform the data
            transformed_list = convert_list_to_mapped_values(value_array, mapping_dict)

            # loop through the table_key_dict and the lookup table. When a match is found add the new mapped key to
            # the existing value of the table_key_dict in a new dict
            converted_dict = {}
            for keys, values in key_values_tables_dict.items():
                for map_key, map_value in mapping_dict.items():
                    if keys == map_key:
                        converted_dict[map_value] = key_values_tables_dict[keys]
            print('NEW_DICTT')
            print(converted_dict)

            return transformed_list, converted_dict

        def wrangle_data_and_attributes_to_be_netcdfified(field_name, field_value):
            """

            """
            print(field_name)
            print(field_value)
            # values to parse into NetCDFVariable attributes list. Once passed they become a netcdf variable attribute.
            # lookup_table is later converted to comments.
            list_of_possible_value = ['long_name', 'units', 'dtype', 'database_field_name', 'lookup_table',
                                      'convert_keys_and_data_to_int8']

            logger.debug('-----------------')
            logger.debug("Field Name: " + str(field_name))
            logger.debug("Field Values: " + str(field_value))
            #attributes_dict = {'description': get_field_description(field_value['database_field_name'])}
            converted_data_array = []
            attributes_dict = {}




            for value in list_of_possible_value:
                logger.debug("Value in list_of_possible_value: " + str(value))

                # if the field value is in the list of accepted values then add to attributes dict
                if field_value.get(value):

                    logger.debug("Processing: " + str(value))
                    try:
                        assert field_value.get('lookup_table')
                        lookup_table_dict = self.get_keys_and_values_table(field_value.get('lookup_table'))
                    except:
                        pass
                    # some key values are already int8 and don't need to be converted. Thus a flag is included in the
                    # field_names
                    if value == 'lookup_table':
                        logger.debug('lookup_table is populated for value: ' + str(value))
                        attributes_dict['comments'] = str(lookup_table_dict)

                    if value == 'convert_keys_and_data_to_int8':
                        logger.debug('convert_keys_and_data_to_int8 is TRUE for value: ' + str(value))
                        # get the transformed data if the data is in string form
                        assert lookup_table_dict
                        logger.debug('converting ' + str(value) + 'string keys to int8....')
                        converted_data_list, converted_key_value_dict = handle_key_value_cases(field_value,
                                                                                               lookup_table_dict)
                        logger.debug('adding converted lookup table as variable attribute...')
                        # this replaces ['comments'] values set in the previous if statement.
                        attributes_dict['comments'] = str(converted_key_value_dict)
                        converted_data_array = np.array(converted_data_list, field_value['dtype'])

                    # for all other values, simply add them to attributes_dict
                    else:
                        attributes_dict[value] = field_value[value]
                # if the value isn't in the list of accepted attributes
                else:
                    logger.debug(str(value) + ' is not set as an accepted attribute.')

            logger.debug('attributes_dict' + str(attributes_dict))
            # print('attributes_dict' + str(attributes_dict))




            # if the data array was converted, return it and the attribute dict.
            if len(converted_data_array) > 0:
                  return converted_data_array, attributes_dict
            # else get the non converted data and return it in an numpy array and the and the attribute dict too
            else:

                #np.min(self.nc_output_dataset.variables['Long']),
                    data_array = np.array(get_data(field_value), dtype=field_value['dtype'])
                #data_array = np.array(get_data(field_value, 10000), dtype=field_value['dtype'])
                # if field_value['dtype'] == 'int8':
                #     data_array_min = np.min(data_array)
                #     print("MIN")
                #     print(data_array_min)
                #     data_array_max = np.max(data_array)
                #     print("MAX")

                    return data_array, attributes_dict

        def get_data(field_name_dict):
            """

            :param field_name_dict:
            :return:
            """
            fill_value = field_name_dict['fill_value']

            # call the sql query and assign results into a python list
            formatted_sql = self.sql_strings_dict_from_yaml['get_data'].format(field_name_dict['database_field_name'],
                                                                               self.survey_id)
            self.cursor.execute(formatted_sql)
            print("NULLS SEARCHING")
            variable_list = []
            for i in self.cursor:
                print("i" + str(i))
                #if the variable has a null value, then swap in the fill_value in its place.
                #The fill_value will be recognised by netcdf and masked out.
                if i == "(None,)":
                    variable_list.append(fill_value)
                else:
                    variable_list.append(i[0])  # getting the first index is required. Otherwise each point is within its own tuple.
                # variable_list.append(
                # i[0])  # getting the first index is required. Otherwise each point is within its own tuple.
            return variable_list
                # return as numpy array, with dtype specified in yaml file.
                #return np.array(variable_list, dtype=field_name_dict['dtype'])


        def convert_list_to_mapped_values(list_to_edit, mapping_dict):

            logger.debug('- - - - - - - - - - - - - - - - - -')
            logger.debug('convert_list_to_mapped_values()')
            logger.debug('list_to_edit: ' + str(list_to_edit))
            logger.debug('mapping_dict: ' + str(mapping_dict))
            transformed_list = []

            for l in list_to_edit:

                for key5, value5, in mapping_dict.items():
                    if l == key5:

                        transformed_list.append(mapping_dict.get(key5))
                    else:
                        pass

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



        yield NetCDFVariable(short_name='ga_gravity_metadata',
                              data=0,
                              dimensions=[],  # Scalar
                              fill_value=None,
                              attributes=generate_ga_metadata_dict(),
                              dtype='int8'  # Byte datatype
                              )



        for field_name, field_value in Grav2NetCDFConverter.settings['field_names'].items():
            print("MA DOG: " + str(field_name) + str(field_value))
            if field_value['dtype'] == 'int8':
                field_value['new'] = int(field_value['fill_value'])
                print("YO YO MUTHA D(*FD")
                print(field_value['new'])
            if field_value['dtype'] == 'int32':
                    print("YO YO MUTHA D(*FD")
                    field_value['new'] = int(field_value['fill_value'])
                    print("PO" + str(field_value['fill_value']))
                    print(field_value['new'])
            if field_value['dtype'] == 'float32':
                field_value['new'] = float(field_value['fill_value'])


            data, attributes = wrangle_data_and_attributes_to_be_netcdfified(field_name, field_value)

            yield NetCDFVariable(short_name=field_value['short_name'],
                                data=data,
                                dimensions=['point'],
                                fill_value= field_value['new'],
                                attributes=attributes
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
    #print(sql_strings_dict)
    # execute sql to return surveys to convert to netcdf
    survey_cursor.execute(sql_strings_dict['sql_get_surveyids'])
    survey_id_list = []

    # tidy the survey id strings
    for survey_row in survey_cursor:
        tidy_sur = re.search('\d+', survey_row[0]).group()
        survey_id_list.append(tidy_sur)

    logger.debug('Survey count =', str(len(survey_id_list)))

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
            #logger.info(g2n.nc_output_dataset.file_format)
            #print(g2n.nc_output_dataset.variables[''])
            print(g2n.nc_output_dataset.variables)
            for data in g2n.nc_output_dataset.variables['Locmeth']:
                print(data)
            for data in g2n.nc_output_dataset.variables['Stattype']:
                print(data)
            print(g2n.nc_output_dataset.variables['Nvalue'])

            del g2n
            break

if __name__ == '__main__':
    main()
