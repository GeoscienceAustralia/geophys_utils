'''
Created on 19 Jul. 2018

@author: Alex Ip
'''

import os
import sys
import logging
import psycopg2
import uuid
from datetime import datetime
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
            
        logger.debug('Connected to database {}:{}/{} as {}'.format(self.postgres_host, 
                                                                       self.postgres_port, 
                                                                       self.postgres_dbname,
                                                                       self.postgres_user))    

    def __del__(self):
        '''
        PostgresDatasetMetadataCache class Destructor
        '''
        if self.db_connection:
            logger.debug('Disconnecting from database {}:{}/{}'.format(self.postgres_host, 
                                                                       self.postgres_port, 
                                                                       self.postgres_dbname))
            self.db_connection.close()
        
        
    def add_dataset(self, dataset):
        '''
        Function to insert or update dataset record
        '''
        # Assign a UUID if one doesn't exist
        if not dataset.metadata_uuid:
            dataset.metadata_uuid = str(uuid.uuid4())
            logger.info('Created new UUID %s' % dataset.metadata_uuid)
            
        #TODO: Do something less hacky and lazy with this
        params = dict(dataset.__dict__)
        
        self.add_survey(ga_survey_id=dataset.ga_survey_id,
                        survey_name=None
                        )
        
        cursor = self.db_connection.cursor()
        
        insert_dataset_sql = '''insert into dataset (
    dataset_title,
    survey_id,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    convex_hull_polygon,
    metadata_uuid
    )
values (
    %(dataset_title)s,
    (select survey_id from survey where ga_survey_id = %(ga_survey_id)s),
    %(longitude_min)s,
    %(longitude_max)s,
    %(latitude_min)s,
    %(latitude_max)s,
    %(convex_hull_polygon)s,
    %(metadata_uuid)s
    )
ON CONFLICT (metadata_uuid) DO NOTHING;
'''
    
        cursor.execute(insert_dataset_sql, params)

        select_dataset_sql = '''select dataset_id
from dataset
where metadata_uuid = %(metadata_uuid)s;
'''
        # Read dataset_id for newly-inserted dataset
        cursor.execute(select_dataset_sql, params)
        dataset_id = next(cursor)[0]
        
        self.add_keywords(dataset_id, dataset.keyword_list)
        self.add_distributions(dataset_id, dataset.distribution_list)


    def add_survey(self, 
                      ga_survey_id,
                      survey_name=None
                      ):
        '''
        Function to insert survey
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
#                  'keywords': keyword_list
                  }
        
        insert_dataset_keyword_sql = """insert into dataset_keyword(dataset_id, keyword_id)
select %(dataset_id)s, 
    keyword_id 
    from keyword
    where keyword_value in ('""" + "', '".join(keyword_list) + """')
ON CONFLICT (dataset_id, keyword_id) DO NOTHING;
"""            
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
        Function to insert new distributions
        '''
        cursor = self.db_connection.cursor()
            
        # Only add each protocol once
        for protocol_value in set([distribution.protocol for distribution in distribution_list]):
            params = {'protocol_value': protocol_value
                      }
            
            insert_protocol_sql = '''insert into protocol(protocol_value)
values(%(protocol_value)s)
ON CONFLICT (protocol_value) DO NOTHING;
'''            
            cursor.execute(insert_protocol_sql, params)
            
            if cursor.rowcount:
                logger.info('New protocol "{}" inserted into protocol table'.format(protocol_value))
            else:
                logger.debug('Protocol "{}" already exists in protocol table'.format(protocol_value))
                
            
        for distribution in distribution_list:
            params = {'dataset_id': dataset_id,
                      'distribution_url': distribution.url,
                      'protocol_value': distribution.protocol
                      }
    
            insert_distribution_sql = '''insert into distribution(dataset_id, distribution_url, protocol_id)
values(%(dataset_id)s, 
    %(distribution_url)s,
    (select protocol_id from protocol where protocol_value = %(protocol_value)s)
    )
ON CONFLICT (distribution_url) DO NOTHING;
'''            
            cursor.execute(insert_distribution_sql, params)
            
            if cursor.rowcount:
                logger.info('Distribution "{}" inserted into table'.format(distribution.url))
            else:
                logger.debug('Distribution "{}" already exists in table'.format(distribution.url))
                

    def search_dataset_distributions(self,
                                     keyword_list,
                                     ll_ur_coords,
                                     protocol
                                     ):
        '''
        Function to return URLs of specified distribution for all datasets with specified keywords and bounding box
        Note that keywords are searched exclusively, i.e. using "and", not "or"
        '''
        cursor = self.db_connection.cursor()
        
        params = {'protocol_value': protocol,
                  'longitude_min': ll_ur_coords[0][0],
                  'longitude_max': ll_ur_coords[1][0],
                  'latitude_min': ll_ur_coords[0][1],
                  'latitude_max': ll_ur_coords[1][1],
                  }

        dataset_search_sql = """select distribution_url
from distribution
inner join protocol using(protocol_id)
inner join dataset using(dataset_id)
"""
        for keyword_index in range(len(keyword_list)):    
            keyword = keyword_list[keyword_index] 
            dataset_search_sql += """inner join (select dataset_id from dataset_keyword
    inner join keyword using(keyword_id)
    where keyword_value = '""" + keyword + """'
    ) keyword{} using(dataset_id)
""".format(keyword_index+1)
        
        dataset_search_sql += """where
    protocol_value = %(protocol_value)s
    and longitude_min <= %(longitude_max)s
    and longitude_max >= %(longitude_min)s
    and latitude_min <= %(latitude_max)s
    and latitude_max >= %(latitude_min)s
"""

        #logger.debug('dataset_search_sql: {}'.format(dataset_search_sql))
        
        cursor.execute(dataset_search_sql, params)
        
        # Return list of distribution_url values
        return [row[0] for row in cursor]
        
def main():
    pdmc = PostgresDatasetMetadataCache(debug=True)
    
    #===========================================================================
    # dataset = Dataset(dataset_title='Test Dataset {}'.format(datetime.now().isoformat()),
    #              ga_survey_id='1',
    #              longitude_min=0,
    #              longitude_max=0,
    #              latitude_min=0,
    #              latitude_max=0,
    #              convex_hull_polygon=None, 
    #              keyword_list=['blah', 'blah blah', 'blah blah blah', 'keyword1'],
    #              distribution_list=[Distribution(url='file://dataset_path',
    #                                              protocol='file'),
    #                                 Distribution(url='https://opendap_endpoint',
    #                                              protocol='opendap'),
    #                                 ],
    #              metadata_uuid=None
    #              )
    # 
    # pdmc.add_dataset(dataset)
    #===========================================================================
    
    print('Search results:')
    for url in pdmc.search_dataset_distributions(keyword_list=['blah', 'blah blah'],
                                                 ll_ur_coords=[[-1, -1], [1, 1]],
                                                 protocol='file'
                                                 ):
        print(url)
                
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
        logger.debug('Logging handlers set up for logger {}'.format(logger.name))

    main()
                
