'''
Created on 12 Dec 2019

@author: Alex Ip
'''
import netCDF4
import logging
import sys
import re
from geophys_utils import NetCDFPointUtils, NetCDFLineUtils, NetCDFGridUtils


logger = logging.getLogger(__name__)


def get_netcdf_util(netcdf_dataset):
    '''
    Function to take a netCDF4 Dataset object, a path to a netCDF file, or an OPeNDAP endpoint
    and return a NetCDFUtils subclass object (i.e. NetCDFPointUtils, NetCDFLineUtils, or NetCDFGridUtils)
    '''
    if type(netcdf_dataset) == str: # String provided as path to netCDF file
        try:
            try:
                _netcdf_dataset = netCDF4.Dataset(netcdf_dataset, 'r')
            except OSError:
                _netcdf_dataset = netCDF4.Dataset(netcdf_dataset + '#fillmismatch', 'r')
            netcdf_dataset = _netcdf_dataset
        except Exception as e:
            logger.error('Unable to open {}: {}'.format(netcdf_dataset, e))
            return 
            
    elif type(netcdf_dataset) != netCDF4.Dataset: # NetCDF4.Dataset object provided
        raise BaseException('Invalid netcdf_dataset type')
    
    # Dataset has line and line_index variables => must be a line dataset
    if set(['line', 'line_index']) <= set(netcdf_dataset.variables.keys()): 
        return NetCDFLineUtils(netcdf_dataset)
    
    # Dataset has 2D (or higher dimensionality) variable with grid_mapping attribute and indexing variables   
    elif len([variable 
              for variable in netcdf_dataset.variables.values()
              if hasattr(variable, 'grid_mapping') # "crs" or similar variable specified
                  and variable.grid_mapping in netcdf_dataset.variables.keys() # "crs" or similar variable exists
                  and len(variable.dimensions) >= 2 # Data variable is a grid
                  and set(variable.dimensions) <= set(netcdf_dataset.variables.keys()) # Indexing variables exist
              ]) > 0:
        return NetCDFGridUtils(netcdf_dataset)
    
    #TODO: Make sure that there are no other tests we could apply here for point datasets
    else:
        return NetCDFPointUtils(netcdf_dataset)
    
    
def main():
    '''
    Quick and dirty test routine. Takes the name of a file containing a list of NCI OPeNDAP endpoints
    and prints each endpoint address along with the class name of the NetCDFUtils subclass
    '''
    nc_list_file_path = sys.argv[1]
    
    with open(nc_list_file_path, 'r') as nc_list_file:
        nc_list = [nc_file_path.strip()
                   for nc_file_path in nc_list_file
                   if re.match('^http.*/dodsC/.*\.nc$', nc_file_path)
                   ]
        nc_list_file.close()
        
        nc_list.sort()
        #logger.debug(nc_list)
        
        for nc in nc_list:
            ncu = get_netcdf_util(nc)
            print(nc, type(ncu).__name__)
            if ncu:
                ncu.close()

if __name__ == '__main__':
    # Setup logging handler if required
    if not logger.handlers:
        # Set handler for root root_logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        #console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

    main()
        