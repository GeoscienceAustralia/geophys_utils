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


class Grav2NetCDFConverter(NetCDFConverter):
    '''
    CSV2NetCDFConverter concrete class for converting CSV data to netCDF
    '''

    field_defs = [
        {'short_name': 'obsno',
         'long_name': 'observation number',
         'column_name': 'obsno',
         'dtype': 'int8',

         },
        {'short_name' : 'grav',
         'long_name' : 'ground gravity',
         'column_name' : 'grav',
         'dtype' : 'float32',
         'units': 'um/s^2' #TODO: Confirm units in DB
        },
        {'short_name': 'latitude',
         'long_name': 'latitude',
         'column_name': 'dlat',
         'dtype': 'float32',
         'units': 'degrees North'
         },
        {'short_name': 'longitude',
         'long_name': 'longitude',
         'column_name': 'dlong',
         'dtype': 'float32',
         'units': 'degrees East'
         },
        {'short_name': 'gndelev',
         'long_name': 'Ground Elevation',
         'column_name': 'GNDELEV',
         'dtype': 'float32',
         'units': 'metres'
         },
        {'short_name': 'meterhgt',
         'long_name': 'Meter Height',
         'column_name': 'meterhgt',
         'dtype': 'float32',
         'units': 'metres'
         },
        {'short_name': 'gridflag',
         'long_name': """Gridding Flag: Flag to note whether data is included in gravity
                        grids produced by GA.
                        Flag = 1 indicated that the data was used in producing GA
                        grids, such as the National Gravity Grids.
                        Flag = 0: the data was not used in GA produced grids. This is
                        usually because the data is less accurate. The data is included
                        as it may still provide useful information.""",
         'column_name': 'gridflag',
         'dtype': 'int8',
         },
        {'short_name': 'reliability',
         'long_name': 'Estimation of Station Reliability',
         'column_name': 'reliab',
         'dtype': 'int8',
         }
    ]



    def __init__(self, nc_out_path, survey_id, con, netcdf_format='NETCDF4_CLASSIC'):
        '''
        Concrete constructor for subclass CSV2NetCDFConverter
        Needs to initialise object with everything that is required for the other Concrete methods
        N.B: Make sure this base class constructor is called from the subclass constructor
        '''

        def get_survey_metadata(survey_id):
            sql_statement = '''select * from gravity.GRAVSURVEYS gs
                         inner join a.surveys using(eno)
                         where gs.surveyid = '201760'
                         and exists (select * from gravity.OBSERVATIONS go
                         where go.surveyid = gs.surveyid
                         and dlong is not null
                         and dlat is not null
                         and status = 'A'
                         and access_code = 'O'
                         and geodetic_datum = 'GDA94'
                         )'''.format(survey_id)
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
                sql_statement = '''select {0} from gravity.OBSERVATIONS go
                                     where gravity.GRAVSURVEYS.surveyid = {1}
                                     and go.dlong is not null
                                     and go.dlat is not null
                                     and status = 'A'
                                     and access_code = 'O'
                                     and geodetic_datum = 'GDA94'
                                     )'''.format(key, survey_id)
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

        # insert survey wide metadata
        metadata_dict = {'title': self.survey_metadata['SURVEYNAME'],
            'Conventions': "CF-1.6,ACDD-1.3",
            'Gravity_Accuracy' : self.survey_metadata['GRAVACC'], #example of how to add oracle fields to global attributes

            'featureType': "trajectory",
            'keywords': 'blah',
            'geospatial_east_min': np.min(self.nc_output_dataset.variables['longitude']),
            'geospatial_east_max': np.max(self.nc_output_dataset.variables['longitude']),
            'geospatial_east_units': "m",
            'geospatial_east_resolution': "point",
            'geospatial_north_min': np.min(self.nc_output_dataset.variables['latitude']),
            'geospatial_north_max': np.min(self.nc_output_dataset.variables['latitude']),
            'geospatial_north_units': "m",
            'geospatial_north_resolution': "point",
            'geospatial_vertical_min': np.min(self.nc_output_dataset.variables[('gndelev')]), # should say if I use gndelev or meter height
            'geospatial_vertical_max': np.max(self.nc_output_dataset.variables[('gndelev')]), # Should this be min(elevation-DOI)?
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
        sql_statement = """select count(*) from gravity.OBSERVATIONS 
where surveyid = '{}'
and dlong is not null
and dlat is not null
and status = 'A'
and access_code = 'O'
and geodetic_datum = 'GDA94' or geodetic_datum = 'WGS84'
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
            and status = 'A'
            and access_code = 'O'
and geodetic_datum = 'GDA94' or geodetic_datum = 'WGS84'
            order by obsno
            """.format(field_def['column_name'], self.survey_id)

            print(sql_statement)
            variable_list = []
            self.cursor.execute(sql_statement)
            for i in self.cursor:
                variable_list.append(
                    i[0])  # getting the first index is required. Otherwise each point is within its own tuple.
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
        gravity_metadata = {key: value.isoformat()
                            if type(value) == datetime
                            else value
                            for key, value in iter(self.survey_metadata.items())
                            if value is not None
                            }
        pprint(gravity_metadata)
        yield NetCDFVariable(short_name='ga_gravity_metadata',
                              data=0,
                              dimensions=[],  # Scalar
                              fill_value=None,
                              attributes=gravity_metadata,
                              dtype='int8'  # Byte datatype
                              )


        for field_def in Grav2NetCDFConverter.field_defs:

            list_of_possible_value = ['long_name', 'units']
            attributes_dict = {}
            for a in list_of_possible_value:
                if field_def.get(a):
                    attributes_dict[a] = field_def[a]
                else:
                    pass


            yield NetCDFVariable(short_name=field_def['short_name'],
                                 data=get_data(field_def),
                                 dimensions=['point'],
                                 fill_value=-1,
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
                        and status = 'A'
                        and access_code = 'O'
                        and geodetic_datum = 'GDA94'
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

        print('Global attributes:')
        pprint(g2n.nc_output_dataset.__dict__)
        print('Dimensions:')
        print(g2n.nc_output_dataset.dimensions)
        print('Variables:')
        print(g2n.nc_output_dataset.variables)
        del g2n
        break

                #count = count + 1


if __name__ == '__main__':
    main()
