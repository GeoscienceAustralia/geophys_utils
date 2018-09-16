'''
Created on 7 Sep. 2018

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
from flask import Flask
from flask_restful import Api
from flask_compress import Compress
from dynamic_kmls import RestfulKMLQuery
from dynamic_kmls.netcdf2kml import settings
import logging

logger = logging.getLogger(__name__)
if settings['global_settings']['debug']:
    logger.setLevel(logging.DEBUG) # Initial logging level for this module
else:
    logger.setLevel(logging.INFO) # Initial logging level for this module

def configure_app_compression(app):
    '''
    Helper function to set app config parameters
    '''
    app.config['COMPRESS_MIMETYPES'] = [RestfulKMLQuery.CONTENT_TYPE,
                                        'text/html', 
                                        'text/css', 
                                        'text/xml', 
                                        'application/json', 
                                        'application/javascript'
                                        ] 
    app.config['COMPRESS_LEVEL'] = 6 
    app.config['COMPRESS_MIN_SIZE'] = 2000 # Don't bother compressing small KML responses 

app = Flask('dynamic_kmls') # Note hard-coded module name

api = Api(app)
api.add_resource(RestfulKMLQuery, '/<string:dataset_type>/query')

if settings['global_settings']['http_compression']:
    configure_app_compression(app)
    Compress(app)

app.run(host='0.0.0.0', debug=settings['global_settings']['debug'])
