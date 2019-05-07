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
#from osgeo.osr import SpatialReference, CoordinateTransformation
import pyproj
import pycrs


# Define CRS name mappings for 
CRS_NAME_MAPPING = {'GDA94': 'EPSG:4283',
                    'EPSG:283': 'EPSG:4283', # EPSG Prefix for UTM zone
                    }

def get_spatial_ref_from_wkt(wkt_or_crs_name):
    '''
    Function to return Proj object for supplied WKT
    @param wkt: Well-known text or CRS name for SpatialReference, including "EPSG:XXXX"
    @return spatial_ref: Proj from WKT
    '''
    #spatial_ref = SpatialReference()
    
    #===========================================================================
    # # Try to resolve WKT
    # result = spatial_ref.ImportFromWkt(wkt_or_crs_name)
    # if not result:
    #     return spatial_ref
    #===========================================================================
    try:
        result = pyproj.Proj(pycrs.parse.from_unknown_text(wkt_or_crs_name).to_proj4())
        return result
    except pycrs.parse.FormatError:
        pass
    
    #===========================================================================
    # # Try to resolve CRS name - either mapped or original
    # result = spatial_ref.SetWellKnownGeogCS(CRS_NAME_MAPPING.get(wkt_or_crs_name) or wkt_or_crs_name) 
    # if not result:
    #     return spatial_ref
    #===========================================================================
    # Try to resolve CRS name - either mapped or original
    try:
        result = pyproj.Proj(pycrs.parse.from_unknown_text(CRS_NAME_MAPPING.get(wkt_or_crs_name) or wkt_or_crs_name).to_proj4())
        return result
    except pycrs.parse.FormatError:
        pass

    # Try common formulations for UTM zones
    #TODO: Fix this so it works in the Northern hemisphere 
    modified_crs_name = re.sub('\s+', '', wkt_or_crs_name.strip().upper())
    utm_match = (re.match('(\w+)/MGAZONE(\d+)', modified_crs_name) or
                 re.match('(\w+)/(\d+)S', modified_crs_name) or
                 re.match('(EPSG:283)(\d{2})', modified_crs_name) 
                 )
    #===========================================================================
    # if utm_match:
    #     modified_crs_name = utm_match.group(1)
    #     utm_zone = int(utm_match.group(2))
    #     result = spatial_ref.SetWellKnownGeogCS(CRS_NAME_MAPPING.get(modified_crs_name) or modified_crs_name)
    # if not result:
    #     spatial_ref.SetUTM(utm_zone, False) # Put this here to avoid potential side effects in downstream code
    #     return spatial_ref
    #===========================================================================
    if utm_match:
        modified_crs_name = utm_match.group(1)
        utm_zone = int(utm_match.group(2))
        try:
            result = pyproj.Proj('+proj=utm +zone={utm_zone} +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'.format(utm_zone=utm_zone))
            return result
        except:
            pass

    assert not result, 'Invalid WKT or CRS name'
    
def get_wkt_from_spatial_ref(spatial_ref):
    '''
    Function to return OGC WKT for supplied Proj instance
    '''
    proj_definition_string =spatial_ref.definition_string()
    
    #TODO: Fix this ugly work-around for missing "+"
    proj_definition_string = re.sub('(\s|^)([^\+])', '\\1+\\2', proj_definition_string)
    
    return pycrs.parse.from_proj4(proj_definition_string).to_ogc_wkt()

def get_coordinate_transformation(from_wkt, to_wkt):
    '''
    Use GDAL to obtain a CoordinateTransformation object to transform coordinates between CRSs or None if no transformation required.
    @parameter from_wkt: WKT or "EPSG:nnnn" string from which to transform
    @parameter to_wkt: WKT or "EPSG:nnnn" string to which to transform
    '''
    # Assume native coordinates if no wkt given
    if from_wkt == to_wkt:
        return None
     
    from_spatial_ref = get_spatial_ref_from_wkt(from_wkt)
    to_spatial_ref = get_spatial_ref_from_wkt(to_wkt)
 
    # This is probably redundant
    if get_wkt_from_spatial_ref(from_spatial_ref) == get_wkt_from_spatial_ref(to_spatial_ref):
        return None
 
    return lambda x, y: pyproj.transform(from_spatial_ref, to_spatial_ref, x, y, z=None, radians=False)


def get_utm_wkt(coordinate, from_wkt):
    '''
    Function to return CRS for UTM zone of specified coordinates.
    Used to transform coords to metres
    @param coordinate: single coordinate pair
    '''
    def get_utm_zone_from_longitude(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def get_hemisphere_from_latitude(latitude):
        if (latitude < 0.0):
            return 'south'
        else:
            return 'north'
        
    latlon_coord_trans = get_coordinate_transformation(
        from_wkt, 'EPSG:4283')
    
    latlon_coord = coordinate if latlon_coord_trans is None else latlon_coord_trans(
        coordinate[0], coordinate[1])
        
    # Set UTM coordinate reference system
    #===========================================================================
    # utm_spatial_ref = SpatialReference()
    # utm_spatial_ref.SetWellKnownGeogCS('WGS84')
    # utm_spatial_ref.SetUTM(utm_getZone(
    #     latlon_coord[0]), utm_isNorthern(latlon_coord[1]))
    #===========================================================================
    utm_spatial_ref = pyproj.Proj('+proj=utm +zone={utm_zone} +{hemisphere} +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'.format(utm_zone=get_utm_zone_from_longitude(latlon_coord[1]),
                                                                                                                                          hemisphere=get_hemisphere_from_latitude(latlon_coord[0])))
    return get_wkt_from_spatial_ref(utm_spatial_ref)

def transform_coords(coordinates, from_wkt, to_wkt):
    '''
    Convert coordinates between specified coordinate reference systems
    @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
    @parameter from_wkt: WKT or "EPSG:nnnn" string from which to transform.
    @parameter to_wkt: WKT or "EPSG:nnnn" string to which to transform.
    '''
    coord_trans = get_coordinate_transformation(
        from_wkt, to_wkt)  # Transform from specified CRS to native CRS

    coordinate_array = np.array(coordinates) # Copy coordinates into fresh array
        
    if not coord_trans:  # No transformation required
        return coordinate_array # Return copy of original coordinates
    
    is_single_coordinate = (coordinate_array.shape == (2,))
    # Reshape 1D array into 2D single coordinate array if only one coordinate provided
    if is_single_coordinate:
        coordinate_array = coordinate_array.reshape((1,2))
        
    new_coordinate_array = np.array(coord_trans(coordinate_array[:,0],
                                                coordinate_array[:,1]))[:,0:2]
    if is_single_coordinate:
        return new_coordinate_array.reshape((2,))
    else: 
        return new_coordinate_array


