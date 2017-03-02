'''
NetCDFUtils class implementing useful functionality against netCDF files

Created on 2Mar.,2017

@author: u76345
'''

import netCDF4
import math
import itertools
import argparse

class NetCDFUtils(object):
    '''
    NetCDFUtils class implementing useful functionality against netCDF files
    '''
    DEFAULT_COPY_OPTIONS = {'complevel': 2, 
                            'zlib': True, 
                            'chunksizes': [8192, 8192]
                            }

    def __init__(self, nc_path):
        self.nc_path = nc_path
        self.netcdf_dataset = netCDF4.Dataset(nc_path, mode="r") 
        
        # Identify all spatial grid variables
        self.data_variable_list = [variable for variable in self.netcdf_dataset.variables.values() 
                                   if hasattr(variable, 'grid_mapping')]
        
        assert self.data_variable_list, 'Unable to determine data variable(s) (must have "grid_mapping" attribute'

                 
    def copy(self, nc_out_path, 
                 datatype_map_dict={},
                 variable_options_dict={},
                 dim_range_dict={},
                 nc_format=None,
                 limit_dim_size=False):
        '''
        Function to copy a netCDF dataset to another one with potential changes to size, format, 
            variable creation options and datatypes.
            
            :param nc_in_path: path to existing netCDF input file 
            :param nc_out_path: path to netCDF output file 
            :param datatype_map_dict: dict containing any maps from source datatype to new datatype.
                e.g. datatype_map_dict={'uint64': 'uint32'}  would onvert all uint64 variables to uint32.
            :param variable_options_dict: dict containing any overrides for per-variable variable creation 
                options. e.g. variable_options_dict={'sst': {'complevel': 2, 'zlib': True}} would apply
                compression to variable 'sst'
            :param dim_range_dict: dict of (start, end+1) tuples keyed by dimension name
            :param nc_format: output netCDF format - 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET', 
                'NETCDF3_64BIT_DATA', 'NETCDF4_CLASSIC', or 'NETCDF4'. Defaults to same as input format.  
            :param limit_dim_size: Boolean flag indicating whether unlimited dimensions should be fixed
        '''       
        
        # Override default variable options with supplied ones for all data variables
        for data_variable in self.data_variable_list:
            variable_dict = dict(NetCDFUtils.DEFAULT_COPY_OPTIONS)
            variable_dict.update(variable_options_dict.get(data_variable.name) or {})
            variable_options_dict[data_variable.name] = variable_dict
                                
        nc_format = nc_format or self.netcdf_dataset.file_format 
        print 'Output format is %s' % nc_format
        
        nc_output_dataset = netCDF4.Dataset(nc_out_path, mode="w", clobber=True, format=nc_format)
        
        try:
            dims_used = set()
            dim_size = {}
            for variable_name, variable in self.netcdf_dataset.variables.iteritems():
                dims_used |= set(variable.dimensions)
                
                for dimension_index in range(len(variable.dimensions)): 
                    dimension_name = variable.dimensions[dimension_index]
                    if dim_size.get(dimension_name):
                        continue
                        
                    dim_range = dim_range_dict.get(dimension_name)
                    if dim_range:
                        dim_size[dimension_name] = dim_range[1] - dim_range[0]
                    else:
                        dim_size[dimension_name] = variable.shape[dimension_index]
                    
            #print dim_size
            
            #Copy dimensions
            for dimension_name, dimension in self.netcdf_dataset.dimensions.iteritems():
                if dimension_name in dims_used: # Discard unused dimensions
                    print 'Copying dimension %s of length %d' % (dimension_name, dim_size[dimension_name])
                    nc_output_dataset.createDimension(dimension_name, 
                                          dim_size[dimension_name] 
                                          if not dimension.isunlimited() or limit_dim_size 
                                          else None)
                else:
                    print 'Skipping unused dimension %s' % dimension_name
    
            # Copy variables
            for variable_name, input_variable in self.netcdf_dataset.variables.iteritems():
                dtype = datatype_map_dict.get(str(input_variable.datatype)) or input_variable.datatype
                
                # Special case for "crs" or "transverse_mercator" - want byte datatype
                if variable_name in ['crs', 'transverse_mercator']: 
                    dtype = 'i1'
                    
                # Start off by copying options from input variable (if specified)
                var_options = input_variable.filters() or {}
                
                # Chunking is defined outside the filters() result
                chunking = input_variable.chunking()
                if chunking and chunking != 'contiguous':
                    # Input variable is chunked
                    input_variable_chunking = [min(chunking[dimension_index], dim_size[input_variable.dimensions[dimension_index]])
                                                 for dimension_index in range(len(chunking))]
                elif len(input_variable.dimensions) == 2: #TODO: Improve this
                    # Input variable is 2D and not chunked - assume row chunking
                    input_variable_chunking = [1, dim_size[input_variable.dimensions[1]]]
                else:
                    # Input variable is not chunked
                    input_variable_chunking = None
                
                # Default to same chunking on input and output
                if input_variable_chunking: 
                    var_options['chunksizes'] = input_variable_chunking
                    
                if hasattr(input_variable, '_FillValue'):
                    var_options['fill_value'] = input_variable._FillValue
                    
                # Apply any supplied options over top of defaults
                var_options.update(variable_options_dict.get(variable_name) or {})
                
                # Ensure chunk sizes aren't bigger than variable sizes
                print var_options
                if var_options.get('chunksizes'):
                    for dimension_index in range(len(input_variable.dimensions)):
                        var_options['chunksizes'][dimension_index] = min(var_options['chunksizes'][dimension_index],
                                                                         dim_size[input_variable.dimensions[dimension_index]])
                             
                options_string = ' with options: %s' % ', '.join(['%s=%s' % item for item in var_options.iteritems()]) if var_options else ''   
                print "Copying variable %s from datatype %s to datatype %s%s" % (variable_name, 
                                                                                 input_variable.datatype, 
                                                                                 dtype, 
                                                                                 options_string)
                # Create output variable
                output_variable = nc_output_dataset.createVariable(variable_name, 
                                              dtype, 
                                              input_variable.dimensions,
                                              **var_options
                                              )
        
                # Copy variable attributes
                print '\tCopying %s attributes: %s' % (variable_name, ', '.join(input_variable.ncattrs()))
                output_variable.setncatts({k: input_variable.getncattr(k) for k in input_variable.ncattrs()})
    
                # Copy data
                if input_variable.shape: # array
                    overall_slices = [slice(*dim_range_dict[input_variable.dimensions[dimension_index]])  
                              if dim_range_dict.get(input_variable.dimensions[dimension_index])
                              else slice(0, input_variable.shape[dimension_index])
                              for dimension_index in range(len(input_variable.dimensions))
                             ]
                    #print overall_slices
                    print '\tCopying %s array data of shape %s' % (variable_name,
                                                                   tuple([overall_slices[dimension_index].stop - overall_slices[dimension_index].start
                                                                          for dimension_index in range(len(input_variable.dimensions))])
                                                                  )
                    
                    if (not input_variable_chunking or 
                        len(input_variable.dimensions) != 2): 
                        # No chunking - Try to copy in one hit
                        output_variable[...] = input_variable[overall_slices]
                    
                    else: # Chunked - perform copy in pieces
                        #TODO: Improve this for small chunks
                        assert len(input_variable.dimensions) == 2, 'Can only chunk copy 2D data at the moment'
                        
                        # Use largest chunk sizes between input and output
                        piece_sizes = [max(var_options['chunksizes'][dimension_index],
                                          input_variable_chunking[dimension_index])
                                      for dimension_index in range(len(input_variable.dimensions))
                                     ]
    
                        piece_index_ranges = [(overall_slices[dimension_index].start / piece_sizes[dimension_index],
                                              int(math.ceil(float(overall_slices[dimension_index].stop) / piece_sizes[dimension_index]))
                                             )
                                             for dimension_index in range(len(input_variable.dimensions))
                                            ]
                                             
                        piece_counts = [int(math.ceil(float(dim_size[input_variable.dimensions[dimension_index]]) / 
                                                      piece_sizes[dimension_index]))
                                        for dimension_index in range(len(input_variable.dimensions)) 
                                        ]
                    
                        print 'Copying %s pieces of size %s cells' % (' x '.join([str(piece_count) for piece_count in piece_counts]),
                                                                      ' x '.join([str(piece_size) for piece_size in piece_sizes])
                                                                      )
                        
                        # Iterate over every piece
                        for piece_indices in itertools.product(*[range(piece_index_ranges[dimension_index][0], 
                                                                      piece_index_ranges[dimension_index][1])
                                                           for dimension_index in range(len(input_variable.dimensions))
                                                          ]
                                                        ):
                            
                                               
                            print '\tCopying piece %s' % (piece_indices,)
                            
                            piece_read_slices = [slice(max(overall_slices[dimension_index].start,
                                                     piece_indices[dimension_index] * piece_sizes[dimension_index]
                                                    ),                                                 
                                                 min(overall_slices[dimension_index].stop,
                                                     (piece_indices[dimension_index] + 1) * piece_sizes[dimension_index]
                                                    )
                                                )
                                           for dimension_index in range(len(input_variable.dimensions))
                                           ]
                            
                            piece_write_slices = [slice(piece_read_slices[dimension_index].start - overall_slices[dimension_index].start,
                                                  piece_read_slices[dimension_index].stop - overall_slices[dimension_index].start,
                                                 )
                                            for dimension_index in range(len(input_variable.dimensions))
                                           ]
                            
                            #print piece_read_slices, piece_write_slices
                            
                            output_variable[piece_write_slices] = input_variable[piece_read_slices]
                    
                else: # scalar variable - simple copy
                    print '\tCopying %s scalar data' % variable_name
                    output_variable = input_variable
                
            # Copy global attributes  
            print "Copying global attributes: %s" % ', '.join(self.netcdf_dataset.__dict__.keys())
            for item, value in self.netcdf_dataset.__dict__.items():
                if type(value) == str:
                    nc_output_dataset.__setattr__(item, value.encode('utf-8'))
                else:
                    nc_output_dataset.__setattr__(item, value)
                    
            print 'Finished copying netCDF dataset %s to %s.' % (self.nc_path, nc_out_path)
        
        finally:
            nc_output_dataset.close()
            
def main():
    '''
    Main function to take command line parameters, perform CSW query and print required output
    '''
    # Define command line arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-c', '--copy', 
                        dest='do_copy', 
                        action='store_const', 
                        const=True, default=False,
                        help='Copy netCDF files')
    parser.add_argument("--chunking", help="comma-separated list of chunk sizes for each dimension",
                        type=str)
    parser.add_argument("input_path")
    parser.add_argument("output_path")
    
    args = parser.parse_args()
    
    if args.do_copy:
        if args.chunking:
            chunking = [int(chunk_size.strip()) for chunk_size in args.chunking.split(',')]
        else:
            chunking = None
            
    ncu = NetCDFUtils(args.input_path)   
    
    ncu.copy(args.output_path, 
             #datatype_map_dict={},
             variable_options_dict={data_variable.name: {'chunksizes': chunking}
                               for data_variable in ncu.data_variable_list
                               } if chunking else {},
             #dim_range_dict={},
             #nc_format=None,
             #limit_dim_size=False
             )
            
        

if __name__ == '__main__':
    main()