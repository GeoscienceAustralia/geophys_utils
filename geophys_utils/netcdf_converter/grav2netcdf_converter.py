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

@author: Alex Ip
'''
from collections import OrderedDict
# from geophys_utils.netcdf_converter.csv2netcdf_converter import CSV2NetCDFConverter
import numpy as np
import cx_Oracle
from geophys_utils.netcdf_converter import NetCDFConverter, NetCDFVariable
import sys
import re
from datetime import datetime
import netCDF4
from pprint import pprint
import yaml
import os


class Grav2NetCDFConverter(NetCDFConverter):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''


    print("HERE")
    try:
        #settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_settings.yml'))
        print(os.path.splitext(__file__)[0] + '_settings.yml')
        settings = yaml.safe_load(open(os.path.splitext(__file__)[0] + '_settings.yml'))
        print(settings)
    except:
        print("boourns")
        settings = {}

    def get_accuracy_method_keys_and_values(self, table):
        print("GET ACC  METH")
        sql_statement = 'select * from gravity.{}'.format(table)
        query_result = self.cursor.execute(sql_statement)
        print(query_result)
        accuracy_method_keys_and_values_dict = {}
        for s in query_result:
            accuracy_method_keys_and_values_dict[s[0]] = s[1]

        print(accuracy_method_keys_and_values_dict)
        return accuracy_method_keys_and_values_dict
        # return dict(zip(field_names, survey_row)

    gravity_metadata_list = [
        'SURVEYID',
        'STATEGROUP',
        'STATIONS',
        'GRAVACC',
        #'GRAVDATUM',
        # as varaible'GNDELEVACC'
        'GNDELEVMETH',
        # should be survey 'GNDELEVDATUM',
        #'RELIAB',
        #survey 'LAYOUT',
    ]

    def __init__(self, nc_out_path, survey_id, con, netcdf_format='NETCDF4'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''

        print('HERE')
        print(type(Grav2NetCDFConverter.settings))
        print("settings")
        print(Grav2NetCDFConverter.settings['field_names'])
        for key, value in Grav2NetCDFConverter.settings['field_names'].items():

                print(key)
                print(type(key))
                print(value)
                print(type(value))
#                print(value['Lat'])



        def get_survey_metadata(survey_id):
            sql_statement = '''
            select * from gravity.GRAVSURVEYS gs
                inner join a.surveys using(eno)
                where gs.surveyid = {0}
                and exists 
                    (select o1.* from gravity.OBSERVATIONS o1
                    left join gravity.OBSERVATIONS o2
                    on o1.surveyid = o2.surveyid
                    and (o1.entrydate > o2.entrydate OR(o1.entrydate = o2.entrydate and o1.obsno > o2.obsno))
                    and o1.geodetic_datum = o2.geodetic_datum
                    and o1.dlat = o2.dlat
                    and o1.dlong = o2.dlong
                    and o1.access_code = o2.access_code
                    and o1.status = o2.status
                        where
                        o1.surveyid = {0}
                        and o1.status = 'A'
                        and o1.access_code = 'O'
                        and o1.dlat is not null
                        and o1.dlong is not null
                        and o1.grav is not null
                        and o1.gndelev is not null
                        and o1.meterhgt is not null
                        and o1.nvalue is not null
                        and o1.ellipsoidhgt is not null
                        and o1.ellipsoidmeterhgt is not null
                        and o1.eno in (select eno from a.surveys where countryid is null or countryid = 'AUS')
                        and o2.obsno is null)'''.format(self.survey_id)

            # sql_statement = '''select * from gravity.GRAVSURVEYS gs
            #              inner join a.surveys using(eno)
            #              where gs.surveyid = {}
            #              and exists
            #              (select * from gravity.OBSERVATIONS go
            #              where go.surveyid = gs.surveyid
            #              and dlong is not null
            #              and dlat is not null
            #              and status = 'A'
            #              and access_code = 'O'
            #              and geodetic_datum = 'GDA94'
            #              )'''.format(survey_id)

            query_result = self.cursor.execute(sql_statement)
            field_names = [field_desc[0] for field_desc in query_result.description]
            survey_row = next(query_result)
            return dict(zip(field_names, survey_row
                            # [str(field) if field else ''
                            #  for field in survey_row
                            #  ]
                            )
                        )

        def get_survey_metadata_in_obs_table(survey_id):
            columns_to_add = {'LOCCACCUOM' : None }
            for key, value in iter(columns_to_add.items()):
                # sql_statement = '''select {0} from gravity.OBSERVATIONS go
                #                      where gravity.GRAVSURVEYS.surveyid = {1}
                #                      and go.dlong is not null
                #                      and go.dlat is not null
                #                      and status = 'A'
                #                      and access_code = 'O'
                #                      and geodetic_datum = 'GDA94'
                #                      )'''.format(key, survey_id)

                query_result = self.cursor.execute(sql_statement)
                value = next(query_result)
            return value

        NetCDFConverter.__init__(self, nc_out_path, netcdf_format)

        self.cursor = con.cursor()
        self.survey_id = survey_id
        self.survey_metadata = get_survey_metadata(survey_id)


    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''

        print("LAT")
        print(self.nc_output_dataset.variables['Lat'])
        print("LONG")
        print(self.nc_output_dataset.variables['Long'])
        # insert survey wide metadata
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
            'geospatial_vertical_min': np.min(self.nc_output_dataset.variables[('Gndelev')]), # should say if I use gndelev or meter height
            'geospatial_vertical_max': np.max(self.nc_output_dataset.variables[('Gndelev')]), # Should this be min(elevation-DOI)?
            'geospatial_vertical_units': "m",
            'geospatial_vertical_resolution': "point",
            'geospatial_vertical_positive': "up",
           # 'history': 'Converted from .dat file {} using definitions file {}'.format(self.aem_dat_path,
                                                                                  #   self.dfn_path),
            'date_created': datetime.now().isoformat()
            }

        return metadata_dict



    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''
#         sql_statement = """select count(*) from gravity.OBSERVATIONS
# where surveyid = '{}'
# and dlong is not null
# and dlat is not null
# and status = 'A'
# and access_code = 'O'
# and geodetic_datum = 'GDA94' or geodetic_datum = 'WGS84'
# """.format(self.survey_id)

        sql_statement = '''
                       select count(*) from gravity.OBSERVATIONS o1
                       left join gravity.OBSERVATIONS o2
                       on 
                           o1.surveyid = o2.surveyid
                           and (o1.entrydate > o2.entrydate
                           OR(o1.entrydate = o2.entrydate and o1.obsno > o2.obsno))
                           and o1.geodetic_datum = o2.geodetic_datum
                           and o1.dlat = o2.dlat
                           and o1.dlong = o2.dlong
                           and o1.access_code = o2.access_code
                           and o1.status = o2.status
                       where 
                           o1.surveyid = {}
                           and o1.status = 'A'
                           and o1.access_code = 'O'
                           and o2.obsno is null
                           and o1.grav is not null
                           and o1.gndelev is not null
                           and o1.meterhgt is not null
                           and o1.nvalue is not null
                           and o1.ellipsoidhgt is not null
                           and o1.ellipsoidmeterhgt is not null
                           and o1.eno in (select
                       eno from a.surveys where countryid is null or countryid = 'AUS')'''.format(self.survey_id)

        print(sql_statement)
        self.cursor.execute(sql_statement)
        point_count = int(next(self.cursor)[0])


        dimensions = OrderedDict()
        dimensions['point'] = point_count  # number of points per survey
        #dimensions['point'] = 1143
        print(dimensions)
        return dimensions

    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects
        '''
        def get_data(field_def):
#             sql_statement = """select {} from gravity.OBSERVATIONS
#             where surveyid = '{}'
#             and dlong is not null
#             and dlat is not null
#             and status = 'A'
#             and access_code = 'O'
# and geodetic_datum = 'GDA94' or geodetic_datum = 'WGS84'
#             order by obsno
#             """.format(field_def['column_name'], self.survey_id)

            sql_statement = '''
                           select o1.{0} from gravity.OBSERVATIONS o1
                           left join gravity.OBSERVATIONS o2
                           on 
                               o1.surveyid = o2.surveyid
                               and (o1.entrydate > o2.entrydate
                               OR(o1.entrydate = o2.entrydate and o1.obsno > o2.obsno))
                               and o1.geodetic_datum = o2.geodetic_datum
                               and o1.dlat = o2.dlat
                               and o1.dlong = o2.dlong
                               and o1.access_code = o2.access_code
                               and o1.status = o2.status
                           where 
                               o1.surveyid = {1}
                               and o1.status = 'A'
                               and o1.access_code = 'O'
                               and o2.obsno is null
                               and o1.grav is not null
                               and o1.gndelev is not null
                               and o1.meterhgt is not null
                               and o1.nvalue is not null
                               and o1.ellipsoidhgt is not null
                               and o1.ellipsoidmeterhgt is not null
                               and o1.eno in (select
                           eno from a.surveys where countryid is null or countryid = 'AUS')'''\
                .format(field_def['database_field_name'], self.survey_id)

            # print(sql_statement)
            variable_list = []
            self.cursor.execute(sql_statement)
            for i in self.cursor:
                variable_list.append(
                    i[0])  # getting the first index is required. Otherwise each point is within its own tuple.

            # # if the field contains characters it must be handled differently
            # if field_def['dtype'] is "char":
            #     print(field_def['dtype'])
            #     #return variable_list
            #     print('variable lsit')
            #     print(variable_list)
            #     str_out = netCDF4.stringtochar(np.array(variable_list, 'S1'))
            #     str_out2 = str_out.flatten('F')
            #     print('str_out2 lsit')
            #     print(str_out2)
            #     print(str_out2.shape)
            #     #b = np.array([list(word) for word in variable_list])
            #
            #     print(str_out.shape)
            #     #return np.array(str_out, dtype=str_out.dtype.str)
            #
            #
            #     return np.array(str_out2, dtype='S1')

            return np.array(variable_list, dtype=field_def['dtype'])
            # brah = np.char.array(variable_list, unicode=True)
            # print(brah)
            # print("ES")
            # return np.char.array(variable_list, unicode=True)
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
        gravity_metadata = {key: value.isoformat()
                            if type(value) == datetime
                            else value
                            for key, value in iter(self.survey_metadata.items())
                            if value is not None
                            #if value in Grav2NetCDFConverter.gravity_metadata_list

                            }

        # add loccaccuom
        # GNDELEVACCUOM
        # METERHGTUNITS
        # GNDELEVACCUOM
        # METERHGTERRUOM
        # GRAVACCUOM - one survey is all null
        # TCERRMETHOD some are all nulls
        # ELLIPSOIDHGTUNITS - always m
        # ELLIPSOIDHGTMETH
        #
        #
        #

        #pprint(gravity_metadata)
        yield NetCDFVariable(short_name='ga_gravity_metadata',
                              data=0,
                              dimensions=[],  # Scalar
                              fill_value=None,
                              attributes=gravity_metadata,
                              dtype='int8'  # Byte datatype
                              )
        list_of_possible_value = ['long_name', 'units', 'dtype', 'comment']

        for key, value in Grav2NetCDFConverter.settings['field_names'].items():

            attributes_dict = {}
            for a in list_of_possible_value:
                print('attribute_dict')
                print(attributes_dict)
                if value.get(a):
                    attributes_dict[a] = value[a]
                else:
                    pass

        # for field_def in Grav2NetCDFConverter.field_names:
        #     list_of_possible_value = ['long_name', 'units', 'dtype', 'comment']
        #     attributes_dict = {}
        #     for a in list_of_possible_value:
        #         if field_def.get(a):
        #             attributes_dict[a] = field_def[a]
        #         else:
        #             pass


                yield NetCDFVariable(short_name=value['short_name'],
                                 data=get_data(value),
                                 dimensions=['point'],
                                 fill_value=None,
                                 attributes=attributes_dict
                                 )



        # yield NetCDFVariable(short_name='test_data',
        #                      data=np.random.random((self.nc_output_dataset.dimensions['lat'].size,
        #                                             self.nc_output_dataset.dimensions['lon'].size)),
        #                      dimensions=['lat', 'lon'],
        #                      fill_value=0.0,
        #                      attributes={'units': 'random crap',
        #                                  'long_name': 'random numbers between 0 and 1'
        #                                  },
        #                      dtype='float32'
        #                      )
        #
        # return


def main():
    # get user input and connect to oracle
    assert len(sys.argv) >= 4, '....'
    nc_out_path = sys.argv[1]
    u_id = sys.argv[2]
    oracle_database = sys.argv[3]
    pw = sys.argv[4]
    con = cx_Oracle.connect(u_id, pw, oracle_database)
    survey_cursor = con.cursor()

    # get a list of all survey ids
    sql_get_surveyids = """select Surveyid from gravity.GRAVSURVEYS gs
                        where exists (select * from gravity.OBSERVATIONS go
                        where go.surveyid = gs.surveyid
                        and dlong is not null
                        and dlat is not null
                        and status = 'A'
                        and access_code = 'O'
                        and geodetic_datum = 'GDA94'
                        )
                        order by gs.SURVEYID"""

    survey_cursor.execute(sql_get_surveyids)
    survey_id_list = []

    # tidy the survey id strings
    for survey_row in survey_cursor:
        tidy_sur = re.search('\d+', survey_row[0]).group()
        survey_id_list.append(tidy_sur)


    #print(survey_id_list)
    print('Survey count =',len(survey_id_list))
    # Loop throught he survey lists to make a netcdf file based off each one.
    for survey in survey_id_list:
        print(survey)
        g2n = Grav2NetCDFConverter(nc_out_path + str(survey) + '.nc', survey, con)
        g2n.get_accuracy_method_keys_and_values('ACCURACYMETHOD')
        g2n.convert2netcdf()
        print('Finished writing netCDF file {}'.format(nc_out_path))

        print('Global attributes:')
        pprint(g2n.nc_output_dataset.__dict__)
        print('Dimensions:')
        print(g2n.nc_output_dataset.dimensions)
        print('Variables:')
        print(g2n.nc_output_dataset.variables)
        print(g2n.nc_output_dataset.file_format)
        print(g2n.nc_output_dataset.variables["Tcunits"][:])
        #g2n.nc_output_dataset.get_accuracy_method_keys_and_values()
        del g2n
        break


if __name__ == '__main__':

    main()
