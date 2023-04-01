'''
Created on 13 Dec 2019

@author: Alex Ip
'''
import logging
import re
import sys

from geophys_utils import get_netcdf_util

logger = logging.getLogger(__name__)


def main():
    '''
    Quick and dirty test routine. Takes the name of a file containing a list of NCI OPeNDAP endpoints
    and prints each endpoint address along with the class name of the NetCDFUtils subclass
    '''
    assert sys.argv, 'Please supply filename containing a list of NCI OPeNDAP endpoints'
    nc_list_file_path = sys.argv[1]

    with open(nc_list_file_path, 'r') as nc_list_file:
        nc_list = [nc_file_path.strip()
                   for nc_file_path in nc_list_file
                   if re.match('^http.*/dodsC/.*\.nc$', nc_file_path)
                   ]
        nc_list_file.close()

        nc_list.sort()
        # logger.debug(nc_list)

        for nc in nc_list:
            ncu = get_netcdf_util(nc)
            print(nc, type(ncu).__name__)
            if ncu:  # If we have a recognised netCDF file (i.e. grid, point or line)

                # TODO: Do some stuff here

                ncu.close()


if __name__ == '__main__':
    # Setup logging handler if required
    if not logger.handlers:
        # Set handler for root root_logger to standard output
        console_handler = logging.StreamHandler(sys.stdout)
        # console_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.DEBUG)
        console_formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

    main()
