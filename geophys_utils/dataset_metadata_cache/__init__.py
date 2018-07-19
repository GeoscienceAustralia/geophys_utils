'''
Created on 19 Jul. 2018

@author: Alex Ip
'''

import abc
import logging
import os
import yaml

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module

settings = yaml.safe_load(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'settings.yml')))


class Distribution(object):
    '''
    Distribution class definition
    '''
    
    def __init__(self,
                 url,
                 protocol
                 ):
        '''
        Distribution class Constructor
        '''
        self.url = url
        self.protocol = protocol

class Dataset(object):
    '''
    Dataset class definition
    '''
    def __init__(self,
                 dataset_title,
                 ga_survey_id,
                 longitude_min,
                 longitude_max,
                 latitude_min,
                 latitude_max,
                 convex_hull_polygon, 
                 keyword_list,
                 distribution_list,
                 metadata_uuid=None
                 ):
        '''
        Dataset class Constructor
        '''
        self.dataset_title = dataset_title
        self.ga_survey_id = ga_survey_id
        self.longitude_min = longitude_min
        self.longitude_max = longitude_max
        self.latitude_min = latitude_min
        self.latitude_max = latitude_max
        self.convex_hull_polygon = convex_hull_polygon 
        self.metadata_uuid = metadata_uuid
        self.keyword_list = keyword_list
        self.distribution_list = distribution_list   



class DatasetMetadataCache(object):
    '''
    DatasetMetadataCache class definition
    '''
    _db_engine = None

    @abc.abstractmethod
    def __init__(self, debug):
        '''
        DatasetMetadataCache class Constructor
        @parameter debug: Boolean flag indicating whether debug output is required
        '''
        # Initialise and set debug property
        self._debug = None
        self.debug = debug

    @abc.abstractmethod
    def insert_or_update_dataset(self, dataset):
        '''
        Function to insert or update dataset record
        '''

    @abc.abstractmethod
    def get_keyword_id(self, 
                       keyword):
        '''
        DatasetMetadataCache class Constructor
        '''

    @abc.abstractmethod
    def get_distribution_id(self, 
                            distribution):
        '''
        DatasetMetadataCache class Constructor
        '''

    @abc.abstractmethod
    def get_protocol_id(self, protocol_value):
        '''
        DatasetMetadataCache class Constructor
        '''
    @property
    def db_engine(self):
        return type(self)._db_engine
    
    @property
    def debug(self):
        return self._debug
    
    @debug.setter
    def debug(self, debug_value):
        if self._debug != debug_value or self._debug is None:
            self._debug = debug_value
            
            if self._debug:
                logger.setLevel(logging.DEBUG)
                logging.getLogger(self.__module__).setLevel(logging.DEBUG)
            else:
                logger.setLevel(logging.INFO)
                logging.getLogger(self.__module__).setLevel(logging.INFO)
                
        logger.debug('Logger {} set to level {}'.format(logger.name, logger.level))
        logging.getLogger(self.__module__).debug('Logger {} set to level {}'.format(self.__module__, logger.level))
        