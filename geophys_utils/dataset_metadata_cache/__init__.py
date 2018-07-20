'''
Created on 20 Jul. 2018

@author: Alex
'''
from ._dataset_metadata_cache import settings, DatasetMetadataCache, Dataset, Distribution
from ._postgres_dataset_metadata_cache import PostgresDatasetMetadataCache
from ._sqlite_dataset_metadata_cache import SQLiteDatasetMetadataCache

def get_dataset_metadata_cache(db_engine='SQLite', *args, **kwargs):  
    '''
    Class factory function to return subclass of DatasetMetadataCache for specified db_engine
    ''' 
    if db_engine == 'SQLite':
        return SQLiteDatasetMetadataCache(*args, **kwargs)
    elif db_engine == 'Postgres':
        return PostgresDatasetMetadataCache(*args, **kwargs)
    else:
        raise BaseException('Unhandled db_engine "{}"'.format(db_engine))