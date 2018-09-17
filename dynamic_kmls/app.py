'''
Created on 7 Sep. 2018

@author: Andrew Turner & Alex Ip, Geoscience Australia
'''
import sys
from flask import Flask
from flask_restful import Api
from flask_compress import Compress
from dynamic_kmls import settings, RestfulKMLQuery
import logging

logger = logging.getLogger()

if settings['global_settings']['debug']:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)
logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))

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

app.run(debug=settings['global_settings']['debug'])

if __name__ == '__main__':
    # Setup logging handlers if required
    if not logger.handlers:
        # Set handler for root logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        #console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
