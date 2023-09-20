'''
Created on 23 Feb 2020

@author: alex
'''
import argparse
import logging
import sys

from geophys_utils._get_netcdf_util import get_netcdf_util

logger = logging.getLogger()
logger.setLevel(logging.INFO)  # Initial logging level for this module


def main():
    '''
    Main function for calling NetCDFUtils.copy function
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--format",
                        help="NetCDF file format (one of 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA')",
                        type=str, default='NETCDF4')
    parser.add_argument("--chunkspec", help="comma-separated list of <dimension_name>/<chunk_size> specifications",
                        type=str)
    parser.add_argument("--complevel", help="Compression level for chunked variables as an integer 0-9. Default is 4",
                        type=int, default=4)
    parser.add_argument('-d', '--debug', action='store_const', const=True, default=False,
                        help='output debug information. Default is no debug info')
    parser.add_argument("input_path")
    parser.add_argument("output_path")

    args = parser.parse_args()

    if args.chunkspec:
        chunk_spec = {dim_name: int(chunk_size)
                      for dim_name, chunk_size in
                      [chunk_spec_string.strip().split('/') for chunk_spec_string in args.chunkspec.split(',')]}
    else:
        chunk_spec = None

    ncu = get_netcdf_util(args.input_path,
                          debug=args.debug
                          )

    ncu.copy(args.output_path,
             # datatype_map_dict={},
             # Compress all chunked variables
             variable_options_dict={variable_name: {'chunksizes': [chunk_spec.get(dimension)
                                                                   for dimension in variable.dimensions
                                                                   ],
                                                    'zlib': bool(args.complevel),
                                                    'complevel': args.complevel
                                                    }
                                    for variable_name, variable in ncu.netcdf_dataset.variables.items()
                                    if (set(variable.dimensions) & set(chunk_spec.keys()))
                                    } if chunk_spec else {},
             # dim_range_dict={'lat': (5,205),'lon': (5,305)},
             # dim_mask_dict={},
             nc_format=args.format,
             # limit_dim_size=False
             )

    logger.debug('Copy complete')


if __name__ == '__main__':
    console_handler = logging.StreamHandler(sys.stdout)
    # console_handler.setLevel(logging.INFO)
    console_handler.setLevel(logging.DEBUG)
    console_formatter = logging.Formatter('%(name)s: %(message)s')
    console_handler.setFormatter(console_formatter)

    if not logger.handlers:
        logger.addHandler(console_handler)
        logger.debug('Logging handlers set up for {}'.format(logger.name))

    main()
