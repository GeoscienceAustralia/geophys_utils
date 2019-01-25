#!/usr/bin/env python

# ===============================================================================
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
# ===============================================================================
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
import gc


def plot_point_dataset(netcdf_point_utils,
                       variable_to_map,
                       utm_bbox=None,
                       plot_title=None,
                       colour_scheme='binary',
                       point_size=1.5,
                       point_step=1,
                       save_path=False,
                       count_dict=None
                       ):
    print(save_path)
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

    try:
        utm_wkt, utm_coords = netcdf_point_utils.utm_coords(netcdf_point_utils.xycoords)
    except:
        pass

    # print(utm_coords)
    # print("HANDLE MASKED COORDS")
    # utm_coords = np.ma.filled(utm_coords.astype(float), np.nan)
    # print(utm_coords)

    utm_zone = get_spatial_ref_from_wkt(utm_wkt).GetUTMZone()  # -ve for Southern Hemisphere
    southern_hemisphere = (utm_zone < 0)
    utm_zone = abs(utm_zone)
    projection = ccrs.UTM(zone=utm_zone, southern_hemisphere=southern_hemisphere)
    # print('utm_zone = {}'.format(utm_zone))
    # #print(utm_coords)

    variable = netcdf_point_utils.netcdf_dataset.variables[variable_to_map]
    print(variable)
    try:
        variable = np.ma.filled(variable.astype(int), 0)
    except:
        pass
    # Set geographic range of plot
    if utm_bbox is None:
        utm_bbox = [
            np.min(utm_coords[:, 0]),
            np.min(utm_coords[:, 1]),
            np.max(utm_coords[:, 0]),
            np.max(utm_coords[:, 1])
        ]
        spatial_mask = np.ones(shape=variable.shape, dtype='Bool')
    else:
        spatial_mask = np.logical_and(np.logical_and((utm_bbox[0] <= utm_coords[:, 0]),
                                                     (utm_coords[:, 0] <= utm_bbox[2])),
                                      np.logical_and((utm_bbox[1] <= utm_coords[:, 1]),
                                                     (utm_coords[:, 1] <= utm_bbox[3]))
                                      )
        utm_coords = utm_coords[spatial_mask]

    print('{} points in UTM bounding box: {}'.format(np.count_nonzero(spatial_mask), utm_bbox))
    # print(utm_coords)

    # compute area
    range_x = utm_bbox[2] - utm_bbox[0]
    range_y = utm_bbox[3] - utm_bbox[1]
    print("range_x in metres: {}".format(range_x))
    print("range_y in metres: {}".format(range_y))
    colour_array = rescale_array(variable[spatial_mask], 0, 1)
    fig_height = range_y * 0.0005
    fig_width = range_x * 0.0005

    font_size = fig_width * 0.1 + 40
    print("fig_width in metres: {}".format(fig_width))
    print("fig_height in metres: {}".format(fig_height))
    fig = plt.figure(figsize=(fig_width, fig_height))

    point_size = fig_height * fig_width * 0.0001 + 0.4
    print(point_size)

    ax = fig.add_subplot(1, 1, 1, projection=projection)
    if plot_title is None:
        if hasattr(netcdf_point_utils.netcdf_dataset, 'title'):
            plot_title = netcdf_point_utils.netcdf_dataset.title + str(point_size)
        else:
            plot_title = netcdf_point_utils.netcdf_dataset.filepath() + ' ' + str(point_size)
    ax.set_title(plot_title, fontsize=font_size)

    # map_image = cimgt.OSM() # https://www.openstreetmap.org/about
    map_image = cimgt.StamenTerrain()  # http://maps.stamen.com/
    # map_image = cimgt.QuadtreeTiles()
    # print(map_image.__dict__)
    ax.add_image(map_image, 10)

    # Compute and set regular tick spacing

    x_increment = pow(10.0, floor(log10(range_x))) / 2
    y_increment = pow(10.0, floor(log10(range_y))) / 2
    x_ticks = np.arange((utm_bbox[0] // x_increment + 1) * x_increment, utm_bbox[2], x_increment)
    y_ticks = np.arange((utm_bbox[1] // y_increment + 1) * y_increment, utm_bbox[3], y_increment)
    plt.xticks(x_ticks, rotation=45)
    plt.yticks(y_ticks)

    # set the x and y axis labels
    plt.xlabel("Eastings (m)", rotation=0, labelpad=20, fontsize=font_size)
    plt.ylabel("Northings (m)", rotation=90, labelpad=20, fontsize=font_size)

    # See link for possible colourmap schemes: https://matplotlib.org/examples/color/colormaps_reference.html
    cm = plt.cm.get_cmap(colour_scheme)

    # build a scatter plot of the specified data, define marker, spatial reference system, and the chosen colour map type
    sc = ax.scatter(utm_coords[::point_step, 0],
                    utm_coords[::point_step, 1],
                    marker='o',
                    c=colour_array[::point_step],
                    s=point_size,
                    alpha=0.9,
                    transform=projection,
                    cmap=cm
                    )
    # set up text box

    measured = count_dict[0]

    if 1 in count_dict:
        interpolated = count_dict[1]
    else:
        interpolated = 0
    if 2 in count_dict:
        extrapolated = count_dict[2]
    else:
        extrapolated = 0

    textstr = 'Measure Points Count: {0} \n' \
              'Interpolated Points Count: {1} \n' \
              'Extrapolated Points Count: {2} \n'.format(measured, interpolated, extrapolated)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.0, 0.0, textstr, transform=ax.transAxes, fontsize=font_size,
            verticalalignment='top', bbox=props)

    # set the colour bar ticks and labels
    try:  # not all variables have units. These will fail on the try and produce the map without tick labels.
        cb = plt.colorbar(sc, ticks=[0, 1])
        cb.ax.set_yticklabels(
            [str(np.min(variable[spatial_mask])), str(np.max(variable[spatial_mask]))])  # vertically oriented colorbar

        cb.set_label("{} {}".format(variable.long_name, variable.units))
    except:
        pass

    if save_path:
        print("IS IT GETTING HERE?")
        saved = False
        try_counter = 0
        while not saved:
            try:
                try_counter = try_counter + 1
                print("attempt: {}".format(try_counter))
                # print(gc.get_objects())
                plt.savefig(save_path, debug=True)
                # print(gc.get_objects())
                plt.clf()
                # print(gc.get_objects())
                plt.close()
                # print(gc.get_objects())
                # cm.close()
                # print(gc.get_objects())
                # sc.close()
                # print(gc.get_objects())
                # fig.clf()
                print('before garbage collection')
                print(gc.get_count())
                saved = True
                print('saved')
                gc.collect()
                print('after garbage collection')
                print(gc.get_count())
            except:
                pass


    else:
        plt.show()
        plt.clf()
        plt.close()
        cm.close()
        sc.close()
        fig.clf()
        gc.collect()


def main():
    '''
    main function for quick and dirty testing
    '''
    # Create NetCDFPointUtils object for specified netCDF dataset
    netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/201780.nc'
    # netcdf_path = 'E:\\Temp\\gravity_point_test\\201780.nc'

    netcdf_dataset = netCDF4.Dataset(netcdf_path)
    npu = NetCDFPointUtils(netcdf_dataset)

    # Plot spatial subset
    plot_point_dataset(npu,
                       'Bouguer',
                       utm_bbox=[630000, 7980000, 680000, 8030000],
                       colour_scheme='gist_heat',
                       point_size=50
                       )


if __name__ == '__main__':
    main()
