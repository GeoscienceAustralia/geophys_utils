'''
Created on 16/11/2016

@author: Alex Ip
'''

class NetCDFLineUtils(object):
    '''
    NetCDFLineUtils class to do various fiddly things with NetCDF geophysics line data files.
    '''

    def __init__(self, netcdf_dataset):
        '''
        NetCDFLineUtils Constructor
        '''
        self.netcdf_dataset = netcdf_dataset
        
        self.crs_variable = netcdf_dataset.variables['crs'] #TODO: Make this more general
        self.crs = self.crs_variable.spatial_ref
        
        
    def get_polygon(self):
        pass
    
    def get_lines(self, bounds, variables):
        pass
    
    def grid_points(self, bounds, grid_resolution, variables, grid_crs=None):
        pass    
    
    def clean_lines(self, lines):
        pass