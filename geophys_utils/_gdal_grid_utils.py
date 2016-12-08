'''
Created on 8Dec.,2016

@author: u76345
'''
import os
import re
import tempfile
import numpy as np
from owslib.wcs import WebCoverageService
from osgeo import gdal
from geophys_utils._crs_utils import transform_coords

def get_gdal_wcs_dataset(wcs_url):
    clean_url = re.match('http[^?]+', wcs_url).group(0)
    temp_xml_path = os.path.join(tempfile.gettempdir(), re.sub('\W', '_', clean_url) + '.xml')
    
    wcs = WebCoverageService(wcs_url, version='1.0.0')
    variable_name = wcs.contents.keys()[0] # Only work with first variable
    
    xml_string = '''<WCS_GDAL>
  <ServiceURL>%s?</ServiceURL>
  <CoverageName>%s</CoverageName>
</WCS_GDAL>''' % (clean_url, variable_name)

    temp_xml_file = open(temp_xml_path, 'w') 
    temp_xml_file.write(xml_string)    
    temp_xml_file.close()                        

    return gdal.Open(temp_xml_path)

def get_gdal_grid_values(gdal_dataset, sample_points, from_crs, band_no=1):
    '''
    Function to return values at a series of points from a GDAL dataset
    '''
    geotransform = gdal_dataset.GetGeoTransform()
    to_crs = gdal_dataset.GetProjection()
    gdal_band = gdal_dataset.GetRasterBand(band_no)
    
    native_sample_points = transform_coords(sample_points, from_crs, to_crs)
    
    values = []
    for point in native_sample_points:
        indices = (int((point[0] - geotransform[0]) / geotransform[1] + 0.5),
                   int((point[1] - geotransform[3]) / geotransform[5] + 0.5))
        value = gdal_band.ReadAsArray(xoff=indices[0], yoff=indices[1], win_xsize=1, win_ysize=1)[0,0]
        #print point, indices, value
        values.append(value)
        
    return np.array(values)