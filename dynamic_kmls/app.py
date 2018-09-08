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

app = Flask('dynamic_kmls')
api = Api(app)
api.add_resource(RestfulKMLQuery, '/<string:dataset_type>/query')

# parser = argparse.ArgumentParser()
# #
# # parser.add_argument("-d", "--dataset_settings", help="Point the flask server to the correct dataset settings in the "
# #                                                      "yaml file. Options include: ground_gravity, aem.", type=str,
# #                     required=False)
#
# args = parser.parse_args()

app.run(debug=DEBUG)
