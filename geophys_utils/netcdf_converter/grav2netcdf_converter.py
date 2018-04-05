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
from geophys_utils.netcdf_converter import sql_statements_dict, NetCDFVariable
import sys




class Grav2NetCDFConverter(sql_statements_dict):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''

    def __init__(self, nc_out_path, survey_id_list, cursor, netcdf_format='NETCDF4_CLASSIC'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''
        sql_statements_dict.__init__(self, nc_out_path, netcdf_format)


        #sql_statements_dict = {
            #"read_grav": "select Surveyid, Grav from gravity.OBSERVATIONS WHERE OBSERVATIONS.Surveyid = %s" % survey_id}
        self.variable_grav_list = []
        count = 1
        for survey_id in survey_id_list:

            #while count < 20:
                self.variable_grav_list = self.read_variable_data(survey_id, 'grav', cursor)
                print self.variable_grav_list
                nc_fid = self.variable_generator()
                #print nc_fid
                # for netcdf_file in self.variable_generator():
                #     print netcdf_file.next()
                #count = count + 1




            #
            # print "------------------------------------"
            # print survey_id
            # print "------------------------------------"
            #
            # # variable_grav = np.array()
            # variable_grav_list = []
            # self.cursor.execute(sql_statements_dict['get_grav'])
            #
            # for i in self.cursor:
            #     variable_grav_list.append(i)
            #     print i
            #yield variable_grav_list  # and all other variables and atributes for that netcdf file/survey
            # get all the info here and yield it - call it

    def read_variable_data(self, survey_id, column_name, cursor):

        sql_statement = "select Surveyid, Grav from gravity.OBSERVATIONS WHERE OBSERVATIONS.Surveyid = %s" % survey_id
        print sql_statement
        variable_grav_list = []
        cursor.execute(sql_statement)
        for i in cursor:
            variable_grav_list.append(i)
            print i
        return variable_grav_list



    def get_global_attributes(self):
        '''
        Concrete method to return dict of global attribute <key>:<value> pairs
        '''
        # insert survey wide metadata
        return {'title': 'test dataset'}

    def get_dimensions(self):
        '''
        Concrete method to return OrderedDict of <dimension_name>:<dimension_size> pairs
        '''
        dimensions = OrderedDict()
        dimensions['point'] = len(self.variable_grav_list)  # number of points per survey

        return dimensions

    def variable_generator(self):
        '''
        Concrete generator to yield NetCDFVariable objects
        '''
        variable_grav_np = np.array(self.variable_grav_list)
        print variable_grav_np.shape
        yield NetCDFVariable(short_name='grav',
                             data=variable_grav_np,
                             dimensions=['point'],
                             fill_value=-1,
                             attributes={'long_name': 'gravity'},
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
    print pw
    print nc_out_path
    print u_id
    print oracle_database

    con = cx_Oracle.connect(u_id, pw, oracle_database)



    cursor = con.cursor()

    # get the list of surveyids
    sql_get_surveyids = "select Surveyid from gravity.GRAVSURVEYS"
    cursor.execute(sql_get_surveyids)
    survey_id_list = []
    for sur in cursor:
        survey_id_list.append(sur)
        print sur

    for survey in survey_id_list:
        g2n = Grav2NetCDFConverter(nc_out_path + survey, survey, cursor)
        g2n.convert2netcdf()
        try:
            print g2n
        except:
            print "nope

if __name__ == '__main__':
    main()
