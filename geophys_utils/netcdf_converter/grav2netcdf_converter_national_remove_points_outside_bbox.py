import sys
from netCDF4 import Dataset
import numpy as np

assert len(sys.argv) == 3, '...'
nc_in_path = sys.argv[1]
nc_out_path = sys.argv[2]

# Input file
ds_in = Dataset('{}'.format(nc_in_path))

# Output file
ds_out = Dataset('{}'.format(nc_out_path), mode='w', format='NETCDF4')

# Copy dimensions
for dim_name, dim in ds_in.dimensions.items():
    ds_out.createDimension(dim_name, len(dim) if not dim.isunlimited() else None)

# Get the indices of the points that lie outside the national bounding box
lon = ds_in.variables['longitude'][:]
lat = ds_in.variables['latitude'][:]
lon_indices_outside_box = np.where((lon < 110.8) | (lon > 156))[0]
lat_indices_outside_box = np.where((lat < -45) | (lat > -9))[0]
indices_outside_box = np.concatenate((lon_indices_outside_box, lat_indices_outside_box))

base_variable_parameters = {'complevel': 4, 'endian': 'little', 'fletcher32': True, 'shuffle': True, 'zlib': True}

for var_name, var in ds_in.variables.items():

    if var_name == 'crs':
        variable_parameters = base_variable_parameters
    else:
        variable_parameters = base_variable_parameters
        variable_parameters['chunksizes'] = var.chunking()
        variable_parameters['fill_value'] = var._FillValue

    out_var = ds_out.createVariable(var_name, var.datatype, var.dimensions, **variable_parameters)

    if 'point' in var.dimensions:
        out_var[:] = np.delete(var[:], indices_outside_box)
    else:
        out_var[:] = var[:]

    out_var.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs() if attr != '_FillValue'})

ds_out.close()

