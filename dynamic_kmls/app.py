'''
Created on 7 Sep. 2018

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
from flask import Flask
from flask_restful import Api
from dynamic_kmls import RestfulKMLQuery
from dynamic_kmls.netcdf2kml import settings
import logging

logger = logging.getLogger(__name__)
if settings['global_settings']['debug']:
    logger.setLevel(logging.DEBUG) # Initial logging level for this module
else:
    logger.setLevel(logging.INFO) # Initial logging level for this module

app = Flask('dynamic_kmls') # Note hard-coded module name
api = Api(app)
api.add_resource(RestfulKMLQuery, '/<string:dataset_type>/query')

app.run(debug=settings['global_settings']['debug'])
