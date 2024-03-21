import sys
import cx_Oracle
from netCDF4 import Dataset
import yaml
import numpy as np
from datetime import datetime
from geophys_utils import points2convex_hull
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

assert len(sys.argv) == 5, '....'
nc_path = sys.argv[1]
oracle_db = sys.argv[2]
u_id = sys.argv[3]
pw = sys.argv[4]
con = cx_Oracle.connect(u_id, pw, oracle_db)
cursor = con.cursor()

ds = Dataset("{}".format(nc_path), mode='a')

yaml_sql_settings = yaml.safe_load(open('grav2netcdf_converter_national_sql_strings.yml'))
sql_strings_dict = yaml_sql_settings['sql_strings_dict']

start_date_sql = sql_strings_dict['get_national_survey_metadata_order_by_startdate']
start_date_query_result = cursor.execute(start_date_sql)
field_names = [field_desc[0] for field_desc in start_date_query_result.description]
start_date_zipped_data = dict(zip(field_names, start_date_query_result.fetchone()))
earliest_start_date = start_date_zipped_data['STARTDATE']

end_date_sql = sql_strings_dict['get_national_survey_metadata_order_by_enddate_desc']
end_date_query_result = cursor.execute(end_date_sql)
field_names = [field_desc[0] for field_desc in end_date_query_result.description]
end_date_zipped_data = dict(zip(field_names, end_date_query_result.fetchone()))
latest_end_date = end_date_zipped_data['ENDDATE']

print(earliest_start_date)
print(latest_end_date)

metadata_dict = {
            'title': "Australian National Ground Gravity Compilation 2023",
            'Conventions': "CF-1.6,ACDD-1.3",
            'keywords': "points, gravity, geophysical survey, Earth sciences, geophysics, geoscientific Information",
            'geospatial_lon_min': np.min(ds.variables['longitude']),
            'geospatial_lon_max': np.max(ds.variables['longitude']),
            'geospatial_lon_units': "degrees_east",
            'geospatial_lon_resolution': "point",
            'geospatial_lat_min': np.min(ds.variables['latitude']),
            'geospatial_lat_max': np.max(ds.variables['latitude']),
            'geospatial_lat_units': "degrees_north",
            'geospatial_lat_resolution': "point",
            'history': "Pulled from point gravity database at Geoscience Australia",
            'summary': "Compilation of ground gravity surveys acquired in Australia. Data acquired from State and "
                       "National Geological Surveys, Academics, and private companies. Data has gone through Quality "
                       "Control processes. Station spacing ranges from 11km to less than 1km. The accuracy of the data "
                       "varies generally with date of acquisition, later data using high precision gravity meters and "
                       "GPS for better accuracies.",
            'location_accuracy_min': np.min(ds.variables['locacc']),
            'location_accuracy_max': np.max(ds.variables['locacc']),
            'location_accuracy_units': "m",
            'elevation_accuracy_min': np.min(ds.variables['gndelevacc']),
            'elevation_accuracy_max': np.max(ds.variables['gndelevacc']),
            'elevation_accuracy_units': "m",
            'gravity_accuracy_min': np.min(ds.variables['gravacc']),
            'gravity_accuracy_max': np.max(ds.variables['gravacc']),
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
    # Compute convex hull and add GML representation to metadata
    coordinates = np.array(list(zip(ds.variables['longitude'][:],
                                    ds.variables['latitude'][:]
                                    )
                                )
                           )
    if len(coordinates) >= 3:
        convex_hull = points2convex_hull(coordinates)
        metadata_dict['geospatial_bounds'] = 'POLYGON((' + ', '.join([' '.join(
            ['%.4f' % ordinate for ordinate in coordinates]) for coordinates in convex_hull]) + '))'
    elif len(coordinates) == 2:  # Two points - make bounding box
        bounding_box = [[min(coordinates[:, 0]), min(coordinates[:, 1])],
                        [max(coordinates[:, 0]), min(coordinates[:, 1])],
                        [max(coordinates[:, 0]), max(coordinates[:, 1])],
                        [min(coordinates[:, 0]), max(coordinates[:, 1])],
                        [min(coordinates[:, 0]), min(coordinates[:, 1])]
                        ]
        metadata_dict['geospatial_bounds'] = 'POLYGON((' + ', '.join([' '.join(
            ['%.4f' % ordinate for ordinate in coordinates]) for coordinates in bounding_box]) + '))'
    elif len(coordinates) == 1:  # Single point
        # TODO: Check whether this is allowable under ACDD
        metadata_dict['geospatial_bounds'] = 'POINT((' + ' '.join(
            ['%.4f' % ordinate for ordinate in coordinates[0]]) + '))'
except:
    logger.warning('Unable to write global attribute "geospatial_bounds"')

print(metadata_dict)

for attribute_name, attribute_value in iter(metadata_dict.items()):
    setattr(ds, attribute_name, attribute_value or '')

cursor.close()
ds.close()
