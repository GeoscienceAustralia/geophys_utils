'''
Created on 20 Jul. 2018

@author: Alex
'''
from geophys_utils.dataset_metadata_cache import SQLiteDatasetMetadataCache

def main(): 
    sdmc = SQLiteDatasetMetadataCache(debug=True) 

    endpoint_list = sdmc.search_dataset_distributions(keyword_list=['AUS', 'ground digital data', 'gravity', 'geophysical survey', 'points'],
                                                 protocol='opendap',
                                                 #ll_ur_coords=[[-179.9, -90.0], [180.0, 90.0]]
                                                      # if search with no paramters reutrns everything.
                                                      #ll lower left, upper right
                                                 ll_ur_coords=[[138.193588256836, -30.5767288208008], [138.480285644531, -30.1188278198242]]
                                                 )
    
    print('Search results:')
    for url in endpoint_list:
        print(url) # You would do your own thing here.

    print('{} endpoints found.'.format(len(endpoint_list)))
if __name__ == '__main__':
    main()

    