'''
Created on 23Nov.,2016

@author: u76345
'''
import numpy as np
import math
from geophys_utils._crs_utils import get_utm_crs, transform_coords


def line_length(line):
    '''
    Function to return length of line
    @param line: iterable containing two two-ordinate iterables, e.g. 2 x 2 array or 2-tuple of 2-tuples
    
    @return length: Distance between start & end points in native units
    '''
    return math.sqrt(math.pow(
        line[1][0] - line[0][0], 2.0) + math.pow(line[1][0] - line[0][0], 2.0))


def point_along_line(line, distance):
    '''
    Function to return a point the specified distance along the line
    @param line: iterable containing two two-ordinate iterables, e.g. 2 x 2 array or 2-tuple of 2-tuples
    @param distance: Distance along line new point should be
    
    @return point: Coordinates of point along line or None if distance > line length
    '''
    length = line_length(line)
    proportion = distance / length

    if proportion < 0 or proportion > 1:
        return None

    return tuple([line[0][dim_index] + proportion * 
                  (line[1][dim_index] - line[0][dim_index]) for dim_index in range(2)])


def utm_coords(coordinate_array, crs):
    '''
    Function to convert coordinates to the appropriate UTM CRS
    @param coordinate_array: Array of shape (n, 2) or iterable containing coordinate pairs
    
    @return crs: WKT for UTM CRS
    @return coordinate_array: Array of shape (n, 2) containing UTM coordinate pairs 
    '''
    native_centre_coords = (np.nanmean(coordinate_array[:,0]), np.nanmean(coordinate_array[:,1]))
    utm_crs = get_utm_crs(native_centre_coords, crs)
    return utm_crs, np.array(transform_coords(coordinate_array, crs, utm_crs))


def coords2distance(coordinate_array):
    '''
    Function to calculate cumulative distance in metres from native (lon/lat) coordinates
    @param coordinate_array: Array of shape (n, 2) or iterable containing coordinate pairs
    
    @return distance_array: Array of shape (n) containing cumulative distances from first coord
    '''
    coord_count = coordinate_array.shape[0]
    distance_array = np.zeros((coord_count,), coordinate_array.dtype)
    cumulative_distance = 0.0
    distance_array[0] = cumulative_distance
    last_point = coordinate_array[0]
    
    for coord_index in range(1, coord_count):
        point = coordinate_array[coord_index]
        distance = math.sqrt(math.pow(point[0] - last_point[0], 2.0) + math.pow(point[1] - last_point[1], 2.0))
        distance = line_length((point, last_point))
        cumulative_distance += distance
        distance_array[coord_index] = cumulative_distance
        last_point = point
        
    return distance_array
    
    
def sample_transect(transect_vertices, crs, sample_metres):
    '''
    Function to return a list of sample points sample_metres apart along lines between transect vertices
    @param transect_vertices: list or array of transect vertex coordinates
    @param crs: coordinate reference system for transect_vertices
    @param sample_metres: distance between sample points in metres
    '''
    transect_vertex_array = np.array(transect_vertices)
    # print 'transect_vertex_array = %s' % transect_vertex_array
    nominal_utm_crs, utm_transect_vertices = utm_coords(transect_vertex_array, crs)
    # print 'nominal_utm_crs = %s' % nominal_utm_crs
    # print 'utm_transect_vertices = %s' % utm_transect_vertices

    sample_points = []
    residual = 0
    for vertex_index in range(len(utm_transect_vertices) - 1):
        utm_line = (utm_transect_vertices[
                    vertex_index], utm_transect_vertices[vertex_index + 1])
        # print 'utm_line = %s' % (utm_line,)
        utm_line_length = line_length(utm_line)
        # print 'utm_line_length = %s' % utm_line_length

        # Skip lines of infinite length
        if utm_line_length == float('inf'):
            continue

        sample_count = (utm_line_length + residual) // sample_metres
        # print 'sample_count = %s' % sample_count
        if not sample_count:
            residual += utm_line_length
            continue

        if residual:  # Use un-sampled distance from last line
            start_point = point_along_line(
                utm_line, sample_metres - residual)
        else:
            start_point = utm_line[0]  # Start at beginning
        # print 'start_point = %s' % (start_point,)

        # Calculate new residual
        residual = (utm_line_length + residual) % sample_metres
        # print 'residual = %s' % residual

        end_point = point_along_line(utm_line, utm_line_length - residual)
        # print 'end_point = %s' % (end_point,)

        try:
            sample_point_array = np.stack([np.linspace(start_point[dim_index], end_point[
                                          dim_index], sample_count + 1) for dim_index in range(2)]).transpose()
            # print 'sample_point_array.shape = %s' %
            # (sample_point_array.shape,)
        except Exception as e:
            print 'Line sampling failed: %s' % e.message
            residual = 0
            continue

        sample_points += list(sample_point_array)

        # Don't double up end point with next start point
        if (not residual) and (vertex_index <
                               len(utm_transect_vertices) - 1):
            sample_points.pop()

    return transform_coords(
        sample_points, nominal_utm_crs, crs), sample_metres
