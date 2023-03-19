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
Functions to obtain edge points and convex hull vertices from a gridded NetCDF dataset

Created on 12Sep.,2016

@author: Alex Ip
'''
import logging
import math

import numpy as np
import shapely.geometry as geometry
from scipy import ndimage
from scipy.spatial import Delaunay
from shapely.ops import cascaded_union, polygonize

from ._array_pieces import array_pieces

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)  # Initial logging level for this module


def get_grid_edge_points(grid_array, dimension_ordinates, nodata_value, max_bytes=None):
    '''
    Function to return a list of coordinates corresponding to pixels on the edge of data-containing areas of the NetCDF dataset
    Data is retrieved in pieces, and pieces around external edge of the full array are padded outwards with no-data to ensure 
    that edges are detected when data runs all the way to the array edge.
    @param grid_array: 2D netCDF4 dataset variable 
    @param dimension_ordinates: tuple of two arrays containing ordinates for grid_array (e.g. lat/lons)
    @param nodata_value: Value in grid array representing null value
    @param max_bytes: Maximum number of bytes to retrieve in each array piece
    @return edge_points: n x 2 array of edge point coordinates
    '''
    assert len(grid_array.shape) == 2, 'grid_array is not 2D'

    edge_point_list = []  # Complete list of edge points (unknown length)
    for piece_array, array_offset in array_pieces(grid_array, max_bytes=max_bytes):
        dimension_subset = [dimension_ordinates[dim_index][array_offset[dim_index]:array_offset[
                                                                                       dim_index] + piece_array.shape[
                                                                                       dim_index]] for dim_index in
                            range(2)]

        # Convert masked array to plain array
        if isinstance(piece_array, np.ma.core.MaskedArray):
            piece_array = piece_array.data

        unpadded_piece_shape = piece_array.shape
        # logger.debug('array_offset = {}, unpadded_piece_shape = {}, piece_array.size = {}'.format(
        #    array_offset, unpadded_piece_shape, piece_array.size))

        # Convert to padded boolean True=data/False=no-data array
        piece_array = np.pad((piece_array != nodata_value), pad_width=1, mode='constant', constant_values=False)

        piece_slices = tuple([slice(*[
            0 if array_offset[dim_index] == 0 else 1,  # Trim padding away if not on lower limit
            piece_array.shape[dim_index] - 1  # Trim off higher end padding
        ])
                              for dim_index in range(2)
                              ])

        # logger.debug(piece_slices)
        # logger.debug('piece_array = {}'.format(piece_array))
        piece_array = piece_array[piece_slices]  # Trim away padding on internal edges if required
        # logger.debug('piece_array = {}'.format(piece_array))

        # Set upper edges to false to get complete edge detection
        piece_array[-1, :] = False
        piece_array[:, -1] = False
        # logger.debug('piece_array = {}'.format(piece_array))

        edge_mask = (ndimage.filters.maximum_filter(piece_array, size=2, mode='constant', cval=False) !=
                     ndimage.filters.minimum_filter(piece_array, size=2, mode='constant', cval=False))

        # edge_mask = filters.sobel(grid_array).astype(np.bool)

        # logger.debug('edge_mask = {}'.format(edge_mask))

        # Detect indices of data/no-data edge points in an n x 2 array
        edge_point_indices = np.transpose(np.array(np.where(edge_mask)))

        # logger.debug('{}\n{}'.format(edge_point_indices.shape, edge_point_indices))

        # Remove index offset introduced by padding lower edges
        for dim_index in range(2):
            if array_offset[dim_index] == 0:
                # logger.debug('removing offset for dimension {}'.format(dim_index))
                edge_point_indices[:, dim_index] = edge_point_indices[:, dim_index] - 1

        # logger.debug('edge_point_indices.shape = {}\nedge_point_indices = {}'.format(edge_point_indices.shape, edge_point_indices))

        # Discard any false edge points detected on inner piece edges
        edge_point_indices = edge_point_indices[np.where(np.all((edge_point_indices < unpadded_piece_shape), axis=1))]

        # logger.debug('edge_point_indices.shape = {}\nedge_point_indices = {}'.format(edge_point_indices.shape, edge_point_indices))

        if edge_point_indices.shape[0]:  # If any edge points found
            piece_edge_points = np.zeros(edge_point_indices.shape, dimension_ordinates[0].dtype)
            # logger.debug('piece_edge_points = {}'.format(piece_edge_points))
            # TODO: Do something more general here to account for YX or XY dimension order - this is for YX only
            pixel_offset = [0.0, 0.0]
            for dim_index in range(2):
                piece_edge_points[:, 1 - dim_index] = dimension_subset[dim_index][edge_point_indices[:, dim_index]]
                pixel_offset[1 - dim_index] = (dimension_subset[dim_index][1] - dimension_subset[dim_index][0]) / 2.0

            # logger.debug('pixel_offset = {}'.format(pixel_offset))
            # logger.debug('piece_edge_points = {}'.format(piece_edge_points))
            piece_edge_points = piece_edge_points - np.array(pixel_offset)  # Apply half pixel offset
            # logger.debug('offset piece_edge_points = {}'.format(piece_edge_points))
            edge_point_list.append(piece_edge_points)
            # logger.debug('{} edge points found for piece'.format(piece_edge_points.shape[0]))
        # else:
        #    logger.debug('No edge points found for piece')

    return np.concatenate(edge_point_list, axis=0)


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
    if dilation != 0:
        convex_hull = convex_hull.buffer(
            dilation, cap_style=2, join_style=2, mitre_limit=tolerance)

    if tolerance != 0:
        convex_hull = convex_hull.simplify(tolerance)

    # Convert polygon to list
    try:
        return [coordinates for coordinates in convex_hull.exterior.coords]
    except AttributeError:
        return [coordinates for coordinates in convex_hull.coords]


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


# =========================================================================
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
# =========================================================================


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
        a = math.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = math.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = math.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
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
