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
Created on 16Nov.,2016

@author: u76345
'''
import re
import numpy as np
from osgeo.osr import SpatialReference, CoordinateTransformation
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

from osgeo import __version__ as osgeo_version

if osgeo_version >= '3.':
    from osgeo.osr import OAMS_TRADITIONAL_GIS_ORDER

# Define CRS name mappings for 
CRS_NAME_MAPPING = {
    'GDA94': 'EPSG:4283',
    'EPSG:283': 'EPSG:4283', # EPSG Prefix for GDA94 UTM zone
    'MGA': 'EPSG:4283', # EPSG Prefix for GDA94 UTM zone
    'WGS84': 'EPSG:4326',
    'AGD84': 'EPSG:4203',
    'AGD66': 'EPSG:4202',
    }

def get_spatial_ref_from_wkt(wkt_or_crs_name):
    '''
    Function to return SpatialReference object for supplied WKT
    @param wkt: Well-known text or CRS name for SpatialReference, including "EPSG:XXXX"
    @return spatial_ref: SpatialReference from WKT
    '''
    if not wkt_or_crs_name:
        return None
    
    spatial_ref = SpatialReference()
    
    result = spatial_ref.SetFromUserInput(wkt_or_crs_name)
    if not result:
        logger.debug('CRS determined using SpatialReference.SetFromUserInput({})'.format(wkt_or_crs_name))
        return spatial_ref
    
    # Try to resolve WKT
    result = spatial_ref.ImportFromWkt(wkt_or_crs_name)
    if not result:
        logger.debug('CRS determined using SpatialReference.ImportFromWkt({})'.format(wkt_or_crs_name))
        return spatial_ref

    # Try to resolve CRS name - either mapped or original
    modified_crs_name = CRS_NAME_MAPPING.get(wkt_or_crs_name) or wkt_or_crs_name
    result = spatial_ref.SetWellKnownGeogCS(modified_crs_name) 
    if not result:
        logger.debug('CRS determined using SpatialReference.SetWellKnownGeogCS({})'.format(modified_crs_name))
        return spatial_ref
    
    match = re.match('EPSG:(\d+)$', wkt_or_crs_name, re.IGNORECASE)
    if match:
        epsg_code = int(match.group(1))
        result = spatial_ref.ImportFromEPSG(epsg_code)
        if not result:
            logger.debug('CRS determined using SpatialReference.ImportFromEPSG({})'.format(epsg_code))
            return spatial_ref
        

    # Try common formulations for UTM zones
    #TODO: Fix this so it works in the Northern hemisphere 
    modified_crs_name = re.sub('\s+', '', wkt_or_crs_name.strip().upper())
    utm_match = (re.match('(\w+)/MGAZONE(\d+)$', modified_crs_name) or
                 re.match('(\w+)/(\d+)S$', modified_crs_name) or
                 re.match('(EPSG:283)(\d{2})$', modified_crs_name) or
                 re.match('(MGA)(\d{2}$)', modified_crs_name) 
                 )
    if utm_match:
        modified_crs_name = utm_match.group(1)
        modified_crs_name = CRS_NAME_MAPPING.get(modified_crs_name) or modified_crs_name
        utm_zone = int(utm_match.group(2))
        result = spatial_ref.SetWellKnownGeogCS(modified_crs_name)
    if not result:
        spatial_ref.SetUTM(utm_zone, False) # Put this here to avoid potential side effects in downstream code
        logger.debug('UTM CRS determined using SpatialReference.SetWellKnownGeogCS({}) (zone{})'.format(modified_crs_name, utm_zone))
        return spatial_ref

    assert not result, 'Invalid WKT or CRS name: "{}"'.format(wkt_or_crs_name)

def get_wkt_from_spatial_ref(spatial_ref):
    '''
    Function to return OGC WKT for supplied Proj instance
    '''
    return spatial_ref.ExportToWkt()


def get_coordinate_transformation(from_wkt, to_wkt):
    '''
    Use GDAL to obtain a CoordinateTransformation object to transform coordinates between CRSs or None if no transformation required.
    @parameter from_wkt: WKT or "EPSG:nnnn" string from which to transform
    @parameter to_wkt: WKT or "EPSG:nnnn" string to which to transform
    '''
    # Assume native coordinates if no wkt given
    if not to_wkt or from_wkt == to_wkt:
        return None
    
    from_spatial_ref = get_spatial_ref_from_wkt(from_wkt)
    to_spatial_ref = get_spatial_ref_from_wkt(to_wkt)

    # This is probably redundant
    if not to_spatial_ref or from_spatial_ref.ExportToWkt() == to_spatial_ref.ExportToWkt():
        return None
    
    # Hack to make sure that traditional x-y coordinate order is always used
    if osgeo_version >= '3.':
        logger.debug('Setting axis mapping strategy to XY  for GDAL 3.X  using OAMS_TRADITIONAL_GIS_ORDER')
        from_spatial_ref.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)
        to_spatial_ref.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER)

    return CoordinateTransformation(from_spatial_ref, to_spatial_ref)

def get_utm_wkt(coordinate, from_wkt):
    '''
    Function to return CRS for UTM zone of specified coordinates.
    Used to transform coords to metres
    @param coordinate: single coordinate pair
    '''
    def utm_getZone(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def utm_isNorthern(latitude):
        if (latitude < 0.0):
            return 0
        else:
            return 1
        
    coordinate_array = np.array(coordinate).reshape((1,2))

    latlon_coord_trans = get_coordinate_transformation(
        from_wkt, 'EPSG:4326')
    latlon_coord = coordinate if latlon_coord_trans is None else latlon_coord_trans.TransformPoints(
        coordinate_array)[0][0:2]
        
    # Set UTM coordinate reference system
    utm_spatial_ref = SpatialReference()
    utm_spatial_ref.SetWellKnownGeogCS('WGS84')
    utm_spatial_ref.SetUTM(utm_getZone(
        latlon_coord[0]), utm_isNorthern(latlon_coord[1]))

    return utm_spatial_ref.ExportToWkt()

def transform_coords(coordinates, from_wkt, to_wkt):
    '''
    Convert coordinates between specified coordinate reference systems
    @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
    @parameter from_wkt: WKT or "EPSG:nnnn" string from which to transform. Defaults to native NetCDF CRS
    @parameter to_wkt: WKT or "EPSG:nnnn" string to which to transform. Defaults to native NetCDF CRS
    '''
    #TODO: Deal with this in a more high-level way
    POINT_CHUNK_SIZE = 1048576 # Number of points to check at any one time to keep memory usage down
            
    coord_trans = get_coordinate_transformation(
        from_wkt, to_wkt)  # Transform from specified CRS to native CRS

    coordinate_array = np.array(coordinates) # Copy coordinates into fresh array
        
    if not coord_trans:  # No transformation required
        return coordinate_array # Return copy of original coordinates
    
    is_single_coordinate = (coordinate_array.shape == (2,))
    # Reshape 1D array into 2D single coordinate array if only one coordinate provided
    if is_single_coordinate:
        coordinate_array = coordinate_array.reshape((1,2))
        
    chunk_start_index = 0
    transformed_coord_chunk_list = []
    while chunk_start_index < len(coordinates):
        chunk_end_index = min(chunk_start_index + POINT_CHUNK_SIZE, len(coordinates))
        logger.debug('Transforming coordinates {} to {} of {}'.format(chunk_start_index, chunk_end_index-1, len(coordinates)))
        transformed_coord_chunk_list.append(np.array(coord_trans.TransformPoints(coordinate_array[slice(chunk_start_index, chunk_end_index)]))[:,0:2])
        chunk_start_index = chunk_end_index
    new_coordinate_array = np.concatenate(transformed_coord_chunk_list)
        
    #new_coordinate_array = np.array(coord_trans.TransformPoints(coordinate_array))[:,0:2] # Old un-chunked coordinate transformation
    
    if is_single_coordinate:
        return new_coordinate_array.reshape((2,))
    else: 
        return new_coordinate_array

def get_reprojected_bounds(bounds, from_wkt, to_wkt):
    '''
    Function to take a bounding box specified in one CRS and return its smallest containing bounding box in a new CRS
    @parameter bounds: bounding box specified as tuple(xmin, ymin, xmax, ymax) in CRS from_wkt
    @parameter from_wkt: WKT for CRS from which to transform bounds
    @parameter to_wkt: WKT for CRS to which to transform bounds
    
    @return reprojected_bounding_box: bounding box specified as tuple(xmin, ymin, xmax, ymax) in CRS to_wkt
    '''
    #if (to_wkt is None) or (from_wkt is None) or (to_wkt == from_wkt):
    if to_wkt == from_wkt:
        return bounds
    
    # Need to look at all four bounding box corners, not just LL & UR
    original_bounding_box =((bounds[0], bounds[1]), (bounds[2], bounds[1]), (bounds[2], bounds[3]), (bounds[0], bounds[3]))
    reprojected_bounding_box = np.array(transform_coords(original_bounding_box, from_wkt, to_wkt))
    
    return [min(reprojected_bounding_box[:,0]), min(reprojected_bounding_box[:,1]), max(reprojected_bounding_box[:,0]), max(reprojected_bounding_box[:,1])]
            
            

