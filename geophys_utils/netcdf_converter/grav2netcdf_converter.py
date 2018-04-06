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

    def __init__(self, nc_out_path, survey_id, con, netcdf_format='NETCDF4_CLASSIC'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''


        self.variables_cursor = con.cursor()
        self.attributes_cursor = con.cursor()

        self.survey_id = survey_id
        #sql_statements_dict = {
            #"read_grav": "select Surveyid, Grav from gravity.OBSERVATIONS WHERE OBSERVATIONS.Surveyid = %s" % survey_id}
        self.variable_grav_list = []
        self.variable_grav_list = self.read_variable_data('grav')
        self.variable_dlat_list = self.read_variable_data('dlat')
        self.variable_dlong_list = self.read_variable_data('dlong')

        #self.variable_generator()
        #self.nc_output_dataset = self.nc_output_dataset



        NetCDFConverter.__init__(self, nc_out_path, netcdf_format)



    def read_variable_data(self, column_name):

        #sql_statement = "select obsno, Surveyid, grav, dlat, dlong from gravity.OBSERVATIONS WHERE OBSERVATIONS.Surveyid = {1}".format(column_name, '197320')

        sql_statement = "select {0} from gravity.OBSERVATIONS " \
                        "where surveyid = {1} " \
                        "and dlong is not null " \
                        "and dlat is not null and" \
                        " status = 'O' order by obsno".format(column_name, '197320')

        print(sql_statement)
        variable_list = []
        self.variables_cursor.execute(sql_statement)
        for i in self.variables_cursor:
            variable_list.append(i[0]) # getting the first index is required. Otherwise wach point is within its own tuple.
        return variable_list


    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''
        # insert survey wide metadata
        attributes_dict = {"GNDELEVACC" : None,
                           "GNDELEVDATUM" : None}
        for key, return_value in attributes_dict.items():
            sql = "select Surveyid, {} from gravity.GRAVSURVEYS where GRAVSURVEYS.Surveyid = {}".format(key, self.survey_id)

            self.attributes_cursor.execute(sql)
            for s in self.attributes_cursor:
                attributes_dict["GNDELEVACC"] = s[1]

        #return {"GNDELEVACC" : attributes_dict["GNDELEVACC"]}
        return {"GNDELEVACC": '5'}
        #return {'title': 'test dataset'}

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''
        dimensions = OrderedDict()
        dimensions['point'] = len(self.variable_grav_list)  # number of points per survey
        print(len(self.variable_grav_list))
        return dimensions

    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects
        '''
        print(self.variable_grav_list)
        variable_grav_np = np.array(self.variable_grav_list, dtype='float32')
        print('shape')
        print(variable_grav_np.flags)
        print(variable_grav_np.shape)
        print(variable_grav_np)
        for i in variable_grav_np:
            print(i)
        variable_dlong_np = np.array(self.variable_dlong_list, dtype='float32')
        variable_dlat_np = np.array(self.variable_dlat_list, dtype='float32')

        print("grav data")


        yield NetCDFVariable(short_name='grav',
                                 data=variable_grav_np,
                                 dimensions=['point'],
                                 fill_value=-1,
                                 attributes={'long_name': 'gravity'},
                                 dtype='int32'
                                 )

        yield NetCDFVariable(short_name='dlong',
                                 data=variable_dlong_np,
                                 dimensions=['point'],
                                 fill_value=-1,
                                 attributes={'long_name': 'longitude'},
                                 dtype='int32'
                                 )
        yield NetCDFVariable(short_name='dlat',
                                 data=variable_dlat_np,
                                 dimensions=['point'],
                                 fill_value=-1,
                                 attributes={'long_name': 'latitude'},
                                 dtype='int32'
                                 )

        # Example of crs variable creation for GDA94
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

    # get the list of surveyids
    sql_get_surveyids = "select Surveyid from gravity.GRAVSURVEYS"
    survey_cursor.execute(sql_get_surveyids)
    survey_id_list = []
    for sur in survey_cursor:

        tidy_sur = re.search('\d+', sur[0]).group()

        survey_id_list.append(tidy_sur)

    count =1
    for survey in survey_id_list:
        while count < 6: #limit the number of surveys for testing purposes


            g2n = Grav2NetCDFConverter(nc_out_path + str(survey) + ".nc", survey, con)
            g2n.convert2netcdf()
            print('Finished writing netCDF file {}'.format(nc_out_path))

            print('Dimensions:')
            print(g2n.nc_output_dataset.dimensions)
            print('Variables:')
            print(g2n.nc_output_dataset.variables)

            try:
                print(g2n)
            except:
                print("nope")
            count = count + 1



if __name__ == '__main__':
    main()
