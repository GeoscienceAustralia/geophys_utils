'''
Created on 16Nov.,2016

@author: u76345
'''
import re
from osgeo.osr import SpatialReference, CoordinateTransformation

def get_spatial_ref_from_crs(crs):
    spatial_ref = SpatialReference()
    # Check for EPSG then Well Known Text
    epsg_match = re.match('^EPSG:(\d+)$', crs)
    if epsg_match:
        spatial_ref.ImportFromEPSG(int(epsg_match.group(1)))
    else:  # Assume valid WKT definition
        spatial_ref.ImportFromWkt(crs)
    return spatial_ref

def get_coordinate_transformation(from_crs, to_crs):
    '''
    Use GDAL to obtain a CoordinateTransformation object to transform coordinates between CRSs or None if no transformation required.
    @parameter from_crs: WKT or "EPSG:nnnn" string from which to transform
    @parameter to_crs: WKT or "EPSG:nnnn" string to which to transform
    '''
    # Assume native coordinates if no crs given
    if from_crs == to_crs:
        return None
    
    from_spatial_ref = get_spatial_ref_from_crs(from_crs)
    to_spatial_ref = get_spatial_ref_from_crs(to_crs)

    # This is probably redundant
    if from_spatial_ref.ExportToWkt() == to_spatial_ref.ExportToWkt():
        return None

    return CoordinateTransformation(from_spatial_ref, to_spatial_ref)

def get_utm_crs(coordinate, from_crs):
    '''
    Function to return CRS for UTM zone of specified coordinates.
    Used to transform coords to metres
    '''
    def utm_getZone(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def utm_isNorthern(latitude):
        if (latitude < 0.0):
            return 0
        else:
            return 1

    latlon_coord_trans = get_coordinate_transformation(
        from_crs, 'EPSG:4326')
    latlon_coord = coordinate if latlon_coord_trans is None else latlon_coord_trans.TransformPoint(
        *coordinate)[0:2]

    # Set UTM coordinate reference system
    utm_spatial_ref = SpatialReference()
    utm_spatial_ref.SetWellKnownGeogCS('WGS84')
    utm_spatial_ref.SetUTM(utm_getZone(
        latlon_coord[0]), utm_isNorthern(latlon_coord[1]))

    return utm_spatial_ref.ExportToWkt()

def transform_coords(coordinates, from_crs, to_crs):
    '''
    Convert coordinates between specified coordinate reference systems
    @parameter coordinates: iterable collection of coordinate pairs or single coordinate pair
    @parameter from_crs: WKT or "EPSG:nnnn" string from which to transform. Defaults to native NetCDF CRS
    @parameter to_crs: WKT or "EPSG:nnnn" string to which to transform. Defaults to native NetCDF CRS
    '''
    coord_trans = get_coordinate_transformation(
        from_crs, to_crs)  # Transform from specified CRS to native CRS

    if not coord_trans:  # No transformation required
        return list(coordinates) # Return copy of original coordinates

    try:  # Multiple coordinates
        return [coordinate[0:2]
                for coordinate in coord_trans.TransformPoints(coordinates)]
    except TypeError:  # Single coordinate
        return coord_trans.TransformPoint(*coordinates)[0:2]


