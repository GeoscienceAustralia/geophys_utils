'''
Created on 19 Jul. 2018

@author: Alex Ip
'''

import abc
import os
import logging
import psycopg2
from geophys_utils.dataset_metadata_cache import settings, DatasetMetadataCache, Dataset, Distribution

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module


class PostgresDatasetMetadataCache(DatasetMetadataCache):
    '''
    PostgresDatasetMetadataCache class definition
    '''
    DEFAULT_POSTGRES_AUTOCOMMIT = True

    _db_engine = 'Postgres'

    def __init__(self,
                 postgres_host=None, 
                 postgres_port=None, 
                 postgres_dbname=None, 
                 postgres_user=None, 
                 postgres_password=None, 
                 autocommit=None, 
                 debug=False
                 ):
        '''
        PostgresDatasetMetadataCache class Constructor
        '''
        super(PostgresDatasetMetadataCache, self).__init__(debug)

        self.postgres_host = postgres_host or settings['POSTGRES_SERVER'] 
        self.postgres_port = postgres_port or settings['POSTGRES_PORT'] 
        self.postgres_dbname = postgres_dbname or settings['POSTGRES_DBNAME'] 
        self.postgres_user = postgres_user or settings['POSTGRES_USER'] 
        self.postgres_password = postgres_password or settings['POSTGRES_PASSWORD'] 
        self.autocommit = autocommit if autocommit is not None else PostgresDatasetMetadataCache.DEFAULT_POSTGRES_AUTOCOMMIT
        
        self.db_connection = psycopg2.connect(host=self.postgres_host, 
                                              port=self.postgres_port, 
                                              dbname=self.postgres_dbname, 
                                              user=self.postgres_user, 
                                              password=self.postgres_password)
        
        if self.autocommit:
            self.db_connection.autocommit = True
            self.db_connection.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
        else:
            self.db_connection.autocommit = False
            self.db_connection.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_READ_COMMITTED)
            
            

    def __del__(self):
        '''
        PostgresDatasetMetadataCache class Destructor
        '''
        if self.db_connection:
            logger.debug('Disconnecting from database {}:{}/{}'.format(self.postgres_host, 
                                                                       self.postgres_port, 
                                                                       self.postgres_dbname))
            self.db_connection.close()
        
        
    def insert_or_update_dataset(self, dataset):
        '''
        Function to insert or update dataset record
        '''
        cursor = self.db_connection.cursor()
        
        sql = ''
        
        cursor.execute(sql)


    def get_survey_id(self, 
                      ga_survey_id,
                      survey_name=None
                      ):
        '''
        Function to return primary key of survey, inserting if necessary
        '''
        cursor = self.db_connection.cursor()
        
        params = {'ga_survey_id': ga_survey_id,
                  'survey_name': survey_name
                  }
        
        select_sql = '''select survey_id, ga_survey_id, survey_name from survey
where ga_survey_id = %(ga_survey_id)s;
'''
        # Check for existing record
        cursor.execute(select_sql, params)
        try:
            row = next(cursor)
            params['survey_id'] = row[0]
            
            # Update details if required
            if survey_name and row[2] != survey_name:
                update_sql = '''update survey
set survey_name = %(survey_name)s
where survey_id = %(survey_id)s;
'''
                cursor.execute(update_sql, params)
                
            return params['survey_id'] 
          
        except StopIteration: # No existing record found       
            insert_sql = '''insert into survey(ga_survey_id, survey_name)
values(%(ga_survey_id)s, %(survey_name)s)
ON CONFLICT (ga_survey_id) DO NOTHING;
'''            
            cursor.execute(insert_sql, params)
        
            if cursor.rowcount:
                logger.info('GA survey ID "{}" inserted into survey table for survey with name "{}"'.format(ga_survey_id, survey_name))
                cursor.execute(select_sql, params)
                row = next(cursor)
                return row[0]            
            else:
                raise BaseException('Unable to insert new survey')
            
            
    def get_keyword_id(self, 
                       keyword):
        '''
        Function to return primary key of keyword, inserting if necessary
        '''
        cursor = self.db_connection.cursor()
        
        params = {'keyword': keyword
                  }
        
        select_sql = '''select * from keyword
where keyword = %(keyword)s;
'''
        # Check for existing record
        cursor.execute(select_sql, params)
        try:
            row = next(cursor)
            return row[0]
          
        except StopIteration: # No existing record found       
            insert_sql = '''insert into keyword(keyword)
values(%(keyword)s)
ON CONFLICT (keyword) DO NOTHING;
'''            
            cursor.execute(insert_sql, params)
        
            if cursor.rowcount:
                logger.info('Keyword "{}" inserted into keyword table'.format(keyword))
                cursor.execute(select_sql, params)
                row = next(cursor)
                return row[0]            
            else:
                raise BaseException('Unable to insert new keyword')
            

    def get_distribution_id(self, 
                            distribution):
        '''
        Function to return primary key of distribution, inserting if necessary
        '''

    def get_protocol_id(self, 
                        protocol_value):
        '''
        Function to return primary key of protocol, inserting if necessary
        '''
        