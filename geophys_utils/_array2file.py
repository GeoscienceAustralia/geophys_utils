'''
Created on 20 Oct. 2017

@author: Alex Ip
'''
from osgeo import gdal

def array2file(data_arrays, projection, geotransform, file_path, file_format, dtype=gdal.GDT_Float32):
    '''
    Function to output array to any(?) GDAL supported raster format and return GDAL dataset
    Formats defined in http://www.gdal.org/formats_list.html
    '''
    data_array_shape = data_arrays[0].shape
    assert [data_array.shape for data_array in data_arrays].count(data_array_shape) == len(data_arrays), 'Data arrays are of different shape'
    print data_array_shape
    driver = gdal.GetDriverByName(file_format)
    gdal_dataset = driver.Create(file_path, 
                                 data_array_shape[1], data_array_shape[0], # Array must be ordered yx
                                 len(data_arrays), 
                                 dtype)
    gdal_dataset.SetGeoTransform(geotransform)
    gdal_dataset.SetProjection(projection)
    
    for band_index in range(len(data_arrays)):
        raster_band = gdal_dataset.GetRasterBand(band_index+1)
        raster_band.WriteArray(data_arrays[band_index])
        
    gdal_dataset.FlushCache()
        
    return gdal_dataset    
