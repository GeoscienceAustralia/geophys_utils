'''
Created on 28 Mar. 2018 for quick-and-dirty testing

@author: u76345
'''
from geophys_utils.netcdf_converter.csv2netcdf_converter import CSV2NetCDFConverter

def main():
    nc_out_path = 'C:\\Temp\\test1.nc'
    c2n = CSV2NetCDFConverter(nc_out_path)
    c2n.convert2netcdf()

if __name__ == '__main__':
    main()
