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

def get_spatial_ref_from_wkt(wkt):
    '''
    Function to return SpatialReference object for supplied WKT
    @param wkt: Well-known text for SpatialReference
    @return spatial_ref: SpatialReference from WKT
    '''
    spatial_ref = SpatialReference()
    # Check for EPSG then Well Known Text
    epsg_match = re.match('^EPSG:(\d+)$', wkt)
    if epsg_match:
        spatial_ref.ImportFromEPSG(int(epsg_match.group(1)))
    else:  # Assume valid WKT definition
        spatial_ref.ImportFromWkt(wkt)
    return spatial_ref

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
    if from_spatial_ref.ExportToWkt() == to_spatial_ref.ExportToWkt():
        return None

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
    coord_trans = get_coordinate_transformation(
        from_wkt, to_wkt)  # Transform from specified CRS to native CRS

    coordinate_array = np.array(coordinates) # Copy coordinates into fresh array
        
    if not coord_trans:  # No transformation required
        return coordinate_array # Return copy of original coordinates
    
    is_single_coordinate = (coordinate_array.shape == (2,))
    # Reshape 1D array into 2D single coordinate array if only one coordinate provided
    if is_single_coordinate:
        coordinate_array = coordinate_array.reshape((1,2))
        
    new_coordinate_array = np.array(coord_trans.TransformPoints(coordinate_array))[:,0:2]
    if is_single_coordinate:
        return new_coordinate_array.reshape((2,))
    else: 
        return new_coordinate_array


