#!/usr/bin/env python

#===============================================================================
#    Copyright 2017 Geoscience Australia
# 
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#===============================================================================
'''
Created on 18Jun.,2018

@author: Andrew Turner & Alex Ip
'''
import netCDF4
from geophys_utils import NetCDFPointUtils
import numpy as np
from geophys_utils import get_spatial_ref_from_wkt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from math import log10, floor, pow

def plot_point_dataset(netcdf_point_utils,
                       variable_to_map,
                       utm_bbox=None,
                       plot_title=None,
                       colour_scheme='binary',
                       point_size=10,
                       point_step=1
                       ):
    '''
    Function to plot data points on a map
    @param netcdf_point_utils: NetCDFPointUtils object wrapping a netCDF dataset
    @param variable_to_map: String specifying variable name for colour map
    @param utm_bbox: UTM Bounding box of form [xmin, ymin, xmax, ymax] or None for all points. Default=None
    @param plot_title: String to prefix before dataset title. Default=None for dataset title or dataset basename
    @param colour_scheme: String specifying colour scheme for data points. Default='binary'
    @param point_size: Point size for data points. Default=10
    @param point_step: Point step between plotted points - used to skip points in dense datasets. Default=1 
    '''
    def rescale_array(input_np_array, new_range_min=0, new_range_max=1):
        old_min = input_np_array.min()
        old_range = input_np_array.max() - old_min
        new_range = new_range_max - new_range_min
    
        scaled_np_array = ((input_np_array - old_min) / old_range * new_range) + new_range_min
    
        return scaled_np_array
    
    if plot_title is None:
        if hasattr(netcdf_point_utils.netcdf_dataset, 'title'):
            plot_title = netcdf_point_utils.netcdf_dataset.title
        else:   
            plot_title = netcdf_point_utils.netcdf_dataset.filepath()

    utm_wkt, utm_coords = netcdf_point_utils.utm_coords(netcdf_point_utils.xycoords)

    utm_zone = get_spatial_ref_from_wkt(utm_wkt).GetUTMZone() # -ve for Southern Hemisphere
    southern_hemisphere = (utm_zone < 0)
    utm_zone = abs(utm_zone)
    projection = ccrs.UTM(zone=utm_zone, southern_hemisphere=southern_hemisphere)
    print('utm_zone = {}'.format(utm_zone))
    #print(utm_coords)

    variable = netcdf_point_utils.netcdf_dataset.variables[variable_to_map]

    # Set geographic range of plot
    if utm_bbox is None:
        utm_bbox = [
            np.min(utm_coords[:,0]),
            np.min(utm_coords[:,1]),
            np.max(utm_coords[:,0]),
            np.max(utm_coords[:,1])
            ]
        spatial_mask = np.ones(shape=variable.shape, dtype='Bool')
    else:
        spatial_mask = np.logical_and(np.logical_and((utm_bbox[0] <= utm_coords[:,0]), 
                                                     (utm_coords[:,0] <= utm_bbox[2])), 
                                      np.logical_and((utm_bbox[1] <= utm_coords[:,1]), 
                                                     (utm_coords[:,1] <= utm_bbox[3]))
                                      )
        utm_coords = utm_coords[spatial_mask]
    
    print('{} points in UTM bounding box: {}'.format(np.count_nonzero(spatial_mask), utm_bbox))
    #print(utm_coords)
    
    colour_array = rescale_array(variable[spatial_mask], 0, 1)

    fig = plt.figure(figsize=(30,30))

    ax = fig.add_subplot(1, 1, 1, projection=projection)

    ax.set_title(plot_title)

    #map_image = cimgt.OSM() # https://www.openstreetmap.org/about
    #map_image = cimgt.StamenTerrain() # http://maps.stamen.com/
    map_image = cimgt.QuadtreeTiles()
    #print(map_image.__dict__)
    ax.add_image(map_image, 10)

    # Compute and set regular tick spacing
    range_x = utm_bbox[2] - utm_bbox[0]
    range_y = utm_bbox[3] - utm_bbox[1]       
    x_increment = pow(10.0, floor(log10(range_x))) / 2
    y_increment = pow(10.0, floor(log10(range_y))) / 2        
    x_ticks = np.arange((utm_bbox[0]//x_increment+1)*x_increment, utm_bbox[2], x_increment)
    y_ticks = np.arange((utm_bbox[1]//y_increment+1)*y_increment, utm_bbox[3], y_increment)
    plt.xticks(x_ticks, rotation=45)
    plt.yticks(y_ticks)

    # set the x and y axis labels
    plt.xlabel("Eastings (m)", rotation=0, labelpad=20)
    plt.ylabel("Northings (m)", rotation=90, labelpad=20)

    # See link for possible colourmap schemes: https://matplotlib.org/examples/color/colormaps_reference.html
    cm = plt.cm.get_cmap(colour_scheme)

    # build a scatter plot of the specified data, define marker, spatial reference system, and the chosen colour map type
    sc = ax.scatter(utm_coords[::point_step,0], 
                    utm_coords[::point_step,1], 
                    marker='o', 
                    c=colour_array[::point_step], 
                    s=point_size, 
                    alpha=0.9, 
                    transform=projection, 
                    cmap=cm
                    )

    # set the colour bar ticks and labels
    cb = plt.colorbar(sc, ticks=[0, 1])
    cb.ax.set_yticklabels([str(np.min(variable[spatial_mask])), str(np.max(variable[spatial_mask]))])  # vertically oriented colorbar
    cb.set_label("{} {}".format(variable.long_name, variable.units))

    plt.show()
    

def main():
    '''
    main function for quick and dirty testing
    '''
    # Create NetCDFPointUtils object for specified netCDF dataset
    netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc'
    #netcdf_path = 'E:\\Temp\\gravity_point_test\\201780.nc'
    
    netcdf_dataset = netCDF4.Dataset(netcdf_path)
    npu = NetCDFPointUtils(netcdf_dataset)

    # Plot spatial subset
    plot_point_dataset(npu, 
                       'Bouguer', 
                       utm_bbox=[630000,7980000,680000,8030000],
                       colour_scheme='gist_heat',
                       point_size=50
                       )

if __name__ == '__main__':
    main()
