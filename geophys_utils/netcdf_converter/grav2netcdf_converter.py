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
import netCDF4


class Grav2NetCDFConverter(NetCDFConverter):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''

    field_defs = [
        {'short_name' : 'grav',
         'long_name' : 'gravity',
         'column_name' : 'grav',
         'dtype' : 'float32',
         'units': 'um/s^2' #TODO: Confirm units in DB

        }
    ]



    def __init__(self, nc_out_path, survey_id, con, netcdf_format='NETCDF4_CLASSIC'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''
        NetCDFConverter.__init__(self, nc_out_path, netcdf_format)

        self.cursor = con.cursor()
        self.survey_id = survey_id
        print(self.__dict__)



    # def read_variable_data(self, column_name):
    #
    #     #sql_statement = "select obsno, Surveyid, grav, dlat, dlong from gravity.OBSERVATIONS WHERE OBSERVATIONS.Surveyid = {1}".format(column_name, '197320')
    #
    #     sql_statement = "select {0} from gravity.OBSERVATIONS " \
    #                     "where surveyid = {1} " \
    #                     "and dlong is not null " \
    #                     "and dlat is not null " \
    #                     "and status = 'O'" \
    #                     "order by obsno".format(column_name, self.survey_id)
    #
    #     print(sql_statement)
    #     variable_list = []
    #     self.variables_cursor.execute(sql_statement)
    #     for i in self.variables_cursor:
    #         variable_list.append(i[0]) # getting the first index is required. Otherwise wach point is within its own tuple.
    #     print("variable_list read from oracle")
    #     print(variable_list)
    #     variable_grav_np = np.array(variable_list, dtype='float32')
    #     return variable_grav_np

    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''
        # insert survey wide metadata
        attributes_dict = {"GNDELEVACC" : None,
                           "GNDELEVDATUM" : None}
        #for key, return_value in attributes_dict.items():
            #sql = "select Surveyid, {} from gravity.GRAVSURVEYS where GRAVSURVEYS.Surveyid = {}".format(key, self.survey_id)
            #
            # self.attributes_cursor.execute(sql)
            # for s in self.attributes_cursor:
            #     attributes_dict["GNDELEVACC"] = s[1]

        #return {"GNDELEVACC" : attributes_dict["GNDELEVACC"]}
        #return {"GNDELEVACC": '5'}
        return {'title': 'test dataset'}

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''
        sql_statement = """select count(*) from gravity.OBSERVATIONS 
where surveyid = '{}'
and dlong is not null
and dlat is not null
and status = 'O'
""".format(self.survey_id)

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
            sql_statement = """select {} from gravity.OBSERVATIONS
            where surveyid = '{}'
            and dlong is not null
            and dlat is not null
            and status = 'O'
            order by obsno
            """.format(field_def['column_name'], self.survey_id)

            print(sql_statement)
            variable_list = []
            self.cursor.execute(sql_statement)
            for i in self.cursor:
                variable_list.append(
                    i[0])  # getting the first index is required. Otherwise wach point is within its own tuple.
            print("variable_list read from oracle")
            print(variable_list)
            return np.array(variable_list, dtype=field_def['dtype'])

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

        for field_def in Grav2NetCDFConverter.field_defs:
           yield NetCDFVariable(short_name=field_def['short_name'],
                                 data=get_data(field_def),
                                 dimensions=['point'],
                                 fill_value=-1,
                                 attributes={'long_name': field_def['long_name'],
                                             'units': field_def['units']
                                             },
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
    assert len(sys.argv) >= 4, '....'
    nc_out_path = sys.argv[1]
    u_id = sys.argv[2]
    oracle_database = sys.argv[3]
    pw = sys.argv[4]
    con = cx_Oracle.connect(u_id, pw, oracle_database)

    survey_cursor = con.cursor()
    #sql_get_surveyids = "select Surveyid from gravity.GRAVSURVEYS"
    # get the list of surveyids

    sql_get_surveyids = """select Surveyid from gravity.GRAVSURVEYS gs
                        where exists (select * from gravity.OBSERVATIONS go
                        where go.surveyid = gs.surveyid
                        and dlong is not null
                        and dlat is not null
                        and go.status = 'O'
                        )
                        order by gs.SURVEYID"""

    survey_cursor.execute(sql_get_surveyids)
    survey_id_list = []
    for survey_row in survey_cursor:
        tidy_sur = re.search('\d+', survey_row[0]).group()
        survey_id_list.append(tidy_sur)
        #(tidy_sur)

    count =1
    print(survey_id_list)
    print('Survey count =',len(survey_id_list))
    for survey in survey_id_list:
        print(survey)
        g2n = Grav2NetCDFConverter(nc_out_path + str(survey) + '.nc', survey, con)
        #print('netcdf dict')
        #print(g2n.nc_output_dataset.__dict__)

        g2n.convert2netcdf()
        print('Finished writing netCDF file {}'.format(nc_out_path))

        print('Dimensions:')
        print(g2n.nc_output_dataset.dimensions)
        print('Variables:')
        print(g2n.nc_output_dataset.variables)
        del g2n

                #count = count + 1


if __name__ == '__main__':
    main()
