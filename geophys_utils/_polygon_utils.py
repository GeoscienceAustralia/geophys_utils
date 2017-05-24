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
Functions to obtain edge points and convex hull vertices from a gridded NetCDF dataset

Created on 12Sep.,2016

@author: Alex Ip
'''
import numpy as np
import math
from scipy import ndimage
import shapely.geometry as geometry
from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Delaunay
from geophys_utils._array_pieces import array_pieces


def get_grid_edge_points(grid_array, dimension_ordinates, nodata_value, max_bytes=None):
    '''
    Function to return a list of coordinates corresponding to pixels on the edge of data-containing areas of the NetCDF dataset
    @param grid_array: 2D netCDF4 dataset variable 
    @param dimension_ordinates: tuple of two arrays containing ordinates for grid_array (e.g. lat/lons)
    @param nodata_value: Value in grid array representing null value
    @param max_bytes: Maximum number of bytes to retrieve in each array piece
    '''
    assert len(grid_array.shape) == 2, 'grid_array is not 2D'

    edge_point_list = []  # Complete list of edge points (unknown length)
    for piece_array, array_offset in array_pieces(grid_array, max_bytes=max_bytes):
        dimension_subset = [dimension_ordinates[dim_index][array_offset[dim_index]:array_offset[
            dim_index] + piece_array.shape[dim_index]] for dim_index in range(2)]

        # Convert masked array to plain array
        if isinstance(piece_array, np.ma.core.MaskedArray):
            piece_array = piece_array.data
# print 'array_offset = %s, piece_array.shape = %s, piece_array.size = %s'
# % (array_offset, piece_array.shape, piece_array.size)

        # Convert to boolean True=data/False=no-data array
        piece_array = (piece_array != nodata_value)

        # Detect edges
        edge_ordinates = np.where(ndimage.filters.maximum_filter(piece_array, size=2) !=
                                  ndimage.filters.minimum_filter(piece_array, size=2))

        if edge_ordinates[0].size: # If any edge points found
            piece_edge_points = np.zeros((edge_ordinates[0].size, 2), dimension_ordinates[0].dtype)
            # TODO: Do something more general here to account for YX or XY
            # dimension order
            piece_edge_points[:, 1] = dimension_subset[0][edge_ordinates[0]]
            piece_edge_points[:, 0] = dimension_subset[1][edge_ordinates[1]]
            edge_point_list += list(piece_edge_points)
#            print '%s edge points found' % piece_edge_points.shape[0]
#        else:
#            print 'No edge points found'

    return edge_point_list


def get_netcdf_edge_points(netcdf_dataset, max_bytes=None):
    '''
    Function to return a list of coordinates corresponding to pixels on the edge of data-containing areas of the NetCDF dataset
    @param netcdf_dataset: netCDF4.Dataset object
    @param max_bytes: Maximum number of bytes to retrieve in each array piece
    '''
    # Find variable with "grid_mapping" attribute - assumed to be 2D data
    # variable
    try:
        data_variable = [variable for variable in netcdf_dataset.variables.values(
        ) if hasattr(variable, 'grid_mapping')][0]
    except:
        raise Exception(
            'Unable to determine data variable (must have "grid_mapping" attribute')
# print 'Variable %s has shape %s' % (data_variable.name,
# data_variable.shape)

    assert len(data_variable.dimensions) == 2, '%s is not 2D' % data_variable.name
    dimension_ordinates = [netcdf_dataset.variables[
        data_variable.dimensions[dim_index]] for dim_index in range(2)]
    nodata_value = data_variable._FillValue

    return get_grid_edge_points(data_variable, dimension_ordinates, nodata_value, max_bytes)


def points2convex_hull(point_list, dilation=0, tolerance=0):
    '''
    Function to return a list of vertex coordinates in the convex hull around data-containing areas of a point list
    @param point_list: Iterable containing coordinates from which to compute convex hull
    @param dilation: distance to dilate convex hull
    @param tolerance: distance tolerance for the simplification of the convex hull
    '''
    convex_hull = geometry.MultiPoint(point_list).convex_hull

    # Offset outward by specified dilation and simplify with specified
    # tolerance
    convex_hull = convex_hull.buffer(
        dilation, cap_style=2, join_style=2, mitre_limit=tolerance).simplify(tolerance)

    # Convert polygon to list
    return [coordinates for coordinates in convex_hull.exterior.coords]


def netcdf2convex_hull(netcdf_dataset, max_bytes=None):
    '''
    Function to return a list of vertex coordinates in the convex hull around data-containing areas of the NetCDF dataset
    @param netcdf_dataset: netCDF4.Dataset object
    @param max_bytes: Maximum number of bytes to retrieve in each array piece
    '''
    # Find variable with "GeoTransform" attribute - assumed to be grid mapping
    # variable
    try:
        grid_mapping_variable = [variable for variable in netcdf_dataset.variables.values(
        ) if hasattr(variable, 'GeoTransform')][0]
    except:
        raise Exception(
            'Unable to determine grid mapping variable (must have "GeoTransform" attribute)')
    GeoTransform = [float(
        number) for number in grid_mapping_variable.GeoTransform.strip().split(' ')]
    avg_pixel_size = (abs(GeoTransform[1]) + abs(GeoTransform[5])) / 2.0

    return points2convex_hull(get_netcdf_edge_points(
        netcdf_dataset, max_bytes), avg_pixel_size, avg_pixel_size)


#=========================================================================
# def netcdf2concave_hull(netcdf_dataset, max_bytes=None):
#     '''
#     Function to return a list of vertex coordinates in the convex hull around data-containing areas of the NetCDF dataset
#     @param netcdf_dataset: netCDF4.Dataset object
#     @param max_bytes: Maximum number of bytes to retrieve in each array piece
#     '''
#     # Find variable with "GeoTransform" attribute - assumed to be grid mapping variable
#     try:
#         grid_mapping_variable = [variable for variable in netcdf_dataset.variables.values() if hasattr(variable, 'GeoTransform')][0]
#     except:
#         raise Exception('Unable to determine grid mapping variable (must have "GeoTransform" attribute')
#     GeoTransform = [float(number) for number in grid_mapping_variable.GeoTransform.strip().split(' ')]
#     avg_pixel_size = (abs(GeoTransform[1]) + abs(GeoTransform[5])) / 2.0
#
#     edge_points = get_edge_points(netcdf_dataset, max_bytes)
#
#     #TODO: Compute alpha value dynamically from point density
#     alpha = 1
#
#     return points2alpha_shape(edge_points, alpha, avg_pixel_size, avg_pixel_size)
#=========================================================================


def points2alpha_shape(points, alpha, dilation=0, tolerance=0):
    """
    Compute the alpha shape (concave hull) of a set of points.
    source: http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the gooeyness of the border. Smaller numbers don't fall inward as much as larger numbers. Too large, and you lose everything!
    @param dilation: distance to dilate convex hull
    @param tolerance: distance tolerance for the simplification of the convex hull
    """
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull

    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add((i, j))
        edge_points.append(coords[[i, j]])

    coords = np.array([point.coords[0]
                       for point in points])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0] - pb[0])**2 + (pa[1] - pb[1])**2)
        b = math.sqrt((pb[0] - pc[0])**2 + (pb[1] - pc[1])**2)
        c = math.sqrt((pc[0] - pa[0])**2 + (pc[1] - pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c) / 2.0
        try:
            # Area of triangle by Heron's formula
            area = math.sqrt(s * (s - a) * (s - b) * (s - c))
            try:
                circum_r = a * b * c / (4.0 * area)
                # Here's the radius filter.
                # print circum_r
                if circum_r < 1.0 / alpha:
                    add_edge(edges, edge_points, coords, ia, ib)
                    add_edge(edges, edge_points, coords, ib, ic)
                    add_edge(edges, edge_points, coords, ic, ia)
            except ZeroDivisionError:
                pass
        except ValueError:
            pass

    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    concave_hull = cascaded_union(triangles)

    # Offset outward by specified dilation and simplify with specified
    # tolerance
    concave_hull = concave_hull.buffer(
        dilation, cap_style=2, join_style=2, mitre_limit=tolerance).simplify(tolerance)

    # Convert polygon to list
    return [coordinates for coordinates in concave_hull.exterior.coords]
