'''
Created on 7 Sep. 2018

@author: Andrew Turner
'''
from flask import Flask
from flask_restful import Api
from dynamic_kmls import settings, DEBUG, RestfulKMLQuery
import logging
from pprint import pprint, pformat

logger = logging.getLogger(__name__)
if DEBUG:
    logger.setLevel(logging.DEBUG) # Initial logging level for this module
else:
    logger.setLevel(logging.INFO) # Initial logging level for this module

app = Flask('dynamic_kmls') # Note hard-coded module name
api = Api(app)
api.add_resource(RestfulKMLQuery, '/<string:dataset_type>/query')

app.run(debug=DEBUG)
