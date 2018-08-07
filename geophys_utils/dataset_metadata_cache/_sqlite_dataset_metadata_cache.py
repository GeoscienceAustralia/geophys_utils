'''
Created on 19 Jul. 2018

@author: Alex Ip
'''

import os
import sys
import re
import logging
import sqlite3
import uuid
from datetime import datetime
from geophys_utils.dataset_metadata_cache import settings, DatasetMetadataCache, Dataset, Distribution

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO) # Initial logging level for this module


class SQLiteDatasetMetadataCache(DatasetMetadataCache):
    '''
    SQLiteDatasetMetadataCache class definition
    '''
    DEFAULT_SQLITE_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'data', 'dataset_metadata_cache.sqlite')
    DDL_SQL_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sqlite_dataset_metadata_cache_ddl.sql')
    
    _db_engine = 'SQLite'

    def __init__(self, sqlite_path=None, debug=False, force_recreate=False):
        '''
        SQLiteIGSNDB Constructor
        @parameter sqlite_path: Path to SQLite database file on disk
        @parameter debug: Boolean flag indicating whether debug output is required
        '''
        super(SQLiteDatasetMetadataCache, self).__init__(debug)
        
        self.sqlite_path = sqlite_path or SQLiteDatasetMetadataCache.DEFAULT_SQLITE_FILE 
            
        self.db_connection = None
        
        # Try to establish and test connection
        try:
            if force_recreate:
                raise BaseException('SQLite Database file recreation forced')
            
            if not os.path.isfile(self.sqlite_path):
                raise BaseException('SQLite Database file {} does not exist'.format(self.sqlite_path))
            
            logger.debug('Connecting to existing database {}'.format(self.sqlite_path))
            self.db_connection = sqlite3.connect(self.sqlite_path) 
            
            cursor = self.db_connection.cursor()
            
            # Perform test query
            cursor.execute('select count(*) from dataset')
            
        except BaseException as e: # Need to re-create DB
            logger.debug('DB re-creation: {}'.format(e))
            # Backup original DB
            if self.db_connection:
                logger.debug('Disconnecting from old database {}'.format(self.sqlite_path))
                self.db_connection.close()
                
            backup_sqlite_path = ('_' + re.sub('\W', '-', datetime.now().isoformat())).join(os.path.splitext(self.sqlite_path))
            logger.info('Backing up old database {} to {}'.format(self.sqlite_path, backup_sqlite_path))
            try:
                os.rename(self.sqlite_path, backup_sqlite_path)
            except FileNotFoundError:
                pass
            
            # Create new database
            logger.info('Creating and connecting to fresh database {}'.format(self.sqlite_path))
            self.db_connection = sqlite3.connect(self.sqlite_path)
            cursor = self.db_connection.cursor()
            
            ddl_sql_file = open(SQLiteDatasetMetadataCache.DDL_SQL_FILE, 'r')
            script_sql = ddl_sql_file.read()
            ddl_sql_file.close()
            
            # Strip comments using non-greedy regex substitutions
            script_sql = re.sub('--.*?$', '', re.sub('/\*.*?\*/', '', script_sql, flags=re.RegexFlag.DOTALL), flags=re.RegexFlag.MULTILINE)
            
            ddl_queries = [ddl_query.strip() + ';\n' for ddl_query in script_sql.split(';') if ddl_query.strip()]
            
            logger.info('Executing DDL script {}'.format(SQLiteDatasetMetadataCache.DDL_SQL_FILE))
            for ddl_query in ddl_queries:
                logger.debug('Executing query:\n{}'.format(ddl_query))
                cursor.execute(ddl_query)
                self.db_connection.commit()
                
            cursor = self.db_connection.cursor()
            
            # Perform test query
            cursor.execute('select count(*) from dataset')
            

        logger.debug('Connected to SQLite database {}'.format(self.sqlite_path))    

    def __del__(self):
        '''
        SQLiteDatasetMetadataCache class Destructor
        '''
        if self.db_connection:
            logger.debug('Disconnecting from database {}'.format(self.sqlite_path))
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
                        survey_name=None, # TODO: Need to fill this in at some stage
                        dataset.start_date,
                        dataset.end_date
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
    metadata_uuid,
    point_count
    )
select
    :dataset_title,
    (select survey_id from survey where ga_survey_id = :ga_survey_id),
    :longitude_min,
    :longitude_max,
    :latitude_min,
    :latitude_max,
    :convex_hull_polygon,
    :metadata_uuid,
    :point_count
where not exists (select dataset_id from dataset where metadata_uuid = :metadata_uuid);
'''
        #logger.debug('insert_dataset_sql: {}'.format(insert_dataset_sql))
        #logger.debug('params: {}'.format(params))
        
        cursor.execute(insert_dataset_sql, params)
        self.db_connection.commit()

        if cursor.rowcount:
            logger.info('New dataset "{}" inserted into table'.format(dataset.dataset_title))
        else:
            logger.debug('Dataset "{}" already exists in table'.format(dataset.dataset_title))
            
        select_dataset_sql = '''select dataset_id
from dataset
where metadata_uuid = :metadata_uuid;
'''
        # Read dataset_id for newly-inserted dataset
        cursor.execute(select_dataset_sql, params)
        dataset_id = next(cursor)[0]
        
        self.add_keywords(dataset_id, dataset.keyword_list)
        self.add_distributions(dataset_id, dataset.distribution_list)


    def add_survey(self, 
                      ga_survey_id,
                      survey_name=None,
                      start_date=None,
                      end_date=None
                      ):
        '''
        Function to insert survey
        '''
        if not ga_survey_id:
            return None
        
        cursor = self.db_connection.cursor()
        
        params = {'ga_survey_id': ga_survey_id,
                  'survey_name': survey_name,
                  'start_date': start_date,
                  'end_date': end_date
                  }
        
        insert_survey_sql = '''insert into survey(ga_survey_id, survey_name, start_date, end_date)
select :ga_survey_id, 
    :survey_name,
    :start_date,
    :end_date
where not exists (select survey_id from survey where ga_survey_id = :ga_survey_id);
'''            
        cursor.execute(insert_survey_sql, params)
        self.db_connection.commit()
    
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
                      'keyword_value': keyword
                      }
            
            insert_keyword_sql = '''insert into keyword(keyword_value)
select :keyword_value
where not exists (select keyword_id from keyword where keyword_value = :keyword_value);
'''            
            cursor.execute(insert_keyword_sql, params)
            self.db_connection.commit()
        
            if cursor.rowcount:
                logger.info('New keyword "{}" inserted into table'.format(keyword))
            else:
                logger.debug('Keyword "{}" already exists in table'.format(keyword))
                
                
        # Insert dataset_keyword records in bulk    
        params = {'dataset_id': dataset_id,
#                  'keywords': keyword_list
                  }
        
        insert_dataset_keyword_sql = """insert into dataset_keyword(dataset_id, keyword_id)
select :dataset_id, 
    keyword_id 
    from keyword
    where keyword_value in ('""" + "', '".join(keyword_list) + """')
        and keyword_id not in (select keyword_id from dataset_keyword where dataset_id = :dataset_id);
"""            
        cursor.execute(insert_dataset_keyword_sql, params)
        self.db_connection.commit()
        

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
# where keyword = :keyword;
# '''
#         # Check for existing record
#         cursor.execute(select_sql, params)
#         try:
#             row = next(cursor)
#             return row[0]
#           
#         except StopIteration: # No existing record found       
#             insert_sql = '''insert or ignore into keyword(keyword)
# values(:keyword)
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
#          self.db_connection.commit()   
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
select :protocol_value
where not exists (select protocol_id from protocol where protocol_value = :protocol_value);
'''            
            cursor.execute(insert_protocol_sql, params)
            self.db_connection.commit()
            
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
select :dataset_id, 
    :distribution_url,
    (select protocol_id from protocol where protocol_value = :protocol_value)
where not exists (select distribution_id from distribution where dataset_id = :dataset_id and distribution_url = :distribution_url);
'''            
            cursor.execute(insert_distribution_sql, params)
            self.db_connection.commit()
            
            if cursor.rowcount:
                logger.info('Distribution "{}" inserted into table'.format(distribution.url))
            else:
                logger.debug('Distribution "{}" already exists in table'.format(distribution.url))
                

    def search_dataset_distributions(self,
                                     keyword_list,
                                     protocol,
                                     ll_ur_coords=None
                                     ):
        '''
        Function to return list of tuples containing metadata for all datasets with specified keywords and bounding box
        Note that keywords are searched exclusively, i.e. using "and", not "or"
        Tuples returned are as follows:
            (ga_survey_id, 
            dataset_title,  
            distribution_url, 
            convex_hull_polygon, 
            longitude_min,
            longitude_max,
            latitude_min,
            latitude_max,
            point_count,
            start_date,
            end_date)    
        '''
        cursor = self.db_connection.cursor()
        
        params = {'protocol_value': protocol}

        if ll_ur_coords:
            params.update({'longitude_min': ll_ur_coords[0][0],
                  'longitude_max': ll_ur_coords[1][0],
                  'latitude_min': ll_ur_coords[0][1],
                  'latitude_max': ll_ur_coords[1][1],
                  })

        dataset_search_sql = """select ga_survey_id,
    dataset_title,
    distribution_url,
    convex_hull_polygon,
    longitude_min,
    longitude_max,
    latitude_min,
    latitude_max,
    point_count,
    start_date,
    end_date    
from distribution
inner join protocol using(protocol_id)
inner join dataset using(dataset_id)
inner join survey using(survey_id)
"""
        for keyword_index in range(len(keyword_list)):    
            keyword = keyword_list[keyword_index] 
            dataset_search_sql += """inner join (select dataset_id from dataset_keyword
    inner join keyword using(keyword_id)
    where keyword_value = '""" + keyword + """'
    ) keyword{} using(dataset_id)
""".format(keyword_index+1)
        
        dataset_search_sql += """where
    protocol_value = :protocol_value
"""
        if ll_ur_coords:
            dataset_search_sql += """    and longitude_min <= :longitude_max
    and longitude_max >= :longitude_min
    and latitude_min <= :latitude_max
    and latitude_max >= :latitude_min
"""
        dataset_search_sql += """order by 1"""

        #logger.debug('dataset_search_sql: {}'.format(dataset_search_sql))
        
        logger.debug('dataset_search_sql: {}'.format(dataset_search_sql))
        logger.debug('params: {}'.format(params))
        cursor.execute(dataset_search_sql, params)
        
        return [tuple(row) for row in cursor]
        
def main():
    sdmc = SQLiteDatasetMetadataCache(debug=True)
    
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
    # sdmc.add_dataset(dataset)
    #===========================================================================
    
    print('Search results:')
    for url in sdmc.search_dataset_distributions(keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
                                                 protocol='opendap',
                                                 ll_ur_coords=None
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
                
