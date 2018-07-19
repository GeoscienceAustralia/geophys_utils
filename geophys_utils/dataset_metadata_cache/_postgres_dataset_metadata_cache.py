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



    def add_survey(self, 
                      ga_survey_id,
                      survey_name=None
                      ):
        '''
        Function to insert survey if necessary
        '''
        if not ga_survey_id:
            return None
        
        cursor = self.db_connection.cursor()
        
        params = {'ga_survey_id': ga_survey_id,
                  'survey_name': survey_name
                  }
        
        insert_survey_sql = '''insert into survey(ga_survey_id, survey_name)
values(%(ga_survey_id)s, %(survey_name)s)
ON CONFLICT (ga_survey_id) DO NOTHING;
'''            
        cursor.execute(insert_survey_sql, params)
    
        if cursor.rowcount:
            logger.info('New survey "{}" inserted into table'.format(ga_survey_id))
        else:
            logger.debug('Survey "{}" already exists in table'.format(ga_survey_id))
            
            
    def add_keywords(self,
                     dataset_id, 
                     keyword_list):
        '''
        Function to return primary key of keyword, inserting if necessary
        '''
        cursor = self.db_connection.cursor()
        
        # Try inserting keywords individually
        for keyword in keyword_list:
            params = {'dataset_id': dataset_id,
                      'keyword': keyword
                      }
            
            insert_keyword_sql = '''insert into keyword(keyword_value)
values(%(keyword)s)
ON CONFLICT (keyword_value) DO NOTHING;
'''            
            cursor.execute(insert_keyword_sql, params)
        
            if cursor.rowcount:
                logger.info('New keyword "{}" inserted into table'.format(keyword))
            else:
                logger.debug('Keyword "{}" already exists in table'.format(keyword))
                
                
        # Insert dataset_keyword records in bulk    
        params = {'dataset_id': dataset_id,
                  'keywords': keyword_list
                  }
        
        insert_dataset_keyword_sql = '''insert into dataset_keyword(dataset_id, keyword_id)
select %(dataset_id)s, 
    keyword_id from keyword
    where keyword_value in %(keywords)s
ON CONFLICT (dataset_id, keyword_id) DO NOTHING;
'''            
        cursor.execute(insert_dataset_keyword_sql, params)
        

#===============================================================================
#     def get_keyword_id(self, 
#                        keyword):
#         '''
#         Function to return primary key of keyword, inserting if necessary
#         '''
#         cursor = self.db_connection.cursor()
#         
#         params = {'keyword': keyword
#                   }
#         
#         select_sql = '''select * from keyword
# where keyword = %(keyword)s;
# '''
#         # Check for existing record
#         cursor.execute(select_sql, params)
#         try:
#             row = next(cursor)
#             return row[0]
#           
#         except StopIteration: # No existing record found       
#             insert_sql = '''insert into keyword(keyword)
# values(%(keyword)s)
# ON CONFLICT (keyword) DO NOTHING;
# '''            
#             cursor.execute(insert_sql, params)
#         
#             if cursor.rowcount:
#                 logger.info('Keyword "{}" inserted into keyword table'.format(keyword))
#                 cursor.execute(select_sql, params)
#                 row = next(cursor)
#                 return row[0]            
#             else:
#                 raise BaseException('Unable to insert new keyword')
#             
#===============================================================================

    def add_distributions(self,
                          dataset_id, 
                          distribution_list):
        '''
        Function to insert new distribution
        '''
        cursor = self.db_connection.cursor()
            
        # Only add each protocol once
        for protocol in set([distribution.protocol for distribution in distribution_list]):
            params = {'protocol_name': protocol
                      }
            
            insert_protocol_sql = '''insert into protocol(protocol_name)
values(%(protocol_name)s)
ON CONFLICT (protocol_name) DO NOTHING;
'''            
            cursor.execute(insert_protocol_sql, params)
            
            if cursor.rowcount:
                logger.info('New protocol "{}" inserted into protocol table'.format(distribution.protocol))
            else:
                logger.debug('Protocol "{}" already exists in protocol table'.format(distribution.protocol))
                
            
        for distribution in distribution_list:
            params = {'dataset_id': dataset_id,
                      'distribution_url': distribution.url,
                      'protocol_name': distribution.protocol
                      }
    
            insert_distribution_sql = '''insert into distribution(dataset_id, distribution_url, protocol_id)
values(%(dataset_id)s, 
    %(distribution_url)s
    (select protocol_id from protocol where protocol_name = %(protocol_name)s)
    )
ON CONFLICT (dataset_id, distribution_url, protocol_id) DO NOTHING;
'''            
            cursor.execute(insert_distribution_sql, params)
            
            if cursor.rowcount:
                logger.info('Distribution "{}" inserted into table'.format(distribution.url))
            else:
                logger.debug('Distribution "{}" already exists in table'.format(distribution.url))
