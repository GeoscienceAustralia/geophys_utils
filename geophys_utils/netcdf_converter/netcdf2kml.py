import simplekml
import netCDF4
import numpy as np
from geophys_utils import NetCDFPointUtils, get_spatial_ref_from_wkt
from geophys_utils import points2convex_hull
from geophys_utils import _polygon_utils


# Create NetCDFPointUtils object for specified netCDF dataset
#netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/195105.nc'
#netcdf_path = 'E:\\Temp\\gravity_point_test\\195256.nc'
#195256
netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_netcdf_4\\201780.nc"
#195105
#201780

netcdf_dataset = netCDF4.Dataset(netcdf_path)
npu = NetCDFPointUtils(netcdf_dataset)
print(npu.point_variables)
print(npu.point_count)


survey_title = str(npu.netcdf_dataset.getncattr('title'))
kml = simplekml.Kml()


#npu.netcdf_dataset.getncattr('title')


#kml.newpoint(name="Kirstenbosch", coords=[(140, -35.988862)])  # lon, lat, optional height
# kml.newpoint(name=npu.netcdf_dataset.variables["Obsno"][0],
#              coords=[
#                  (npu.netcdf_dataset.variables["longitude"][0], npu.netcdf_dataset.variables["latitude"][0])])



# Set variables
west_extent = npu.netcdf_dataset.getncattr('geospatial_lon_min')
print("west")
print(west_extent)

east_extent = npu.netcdf_dataset.getncattr('geospatial_lon_max')
print("east")
print(east_extent)


south_extent = npu.netcdf_dataset.getncattr('geospatial_lat_min')
print("south")
print(south_extent)

north_extent = npu.netcdf_dataset.getncattr('geospatial_lat_max')
print("north")
print(north_extent)

#set lod variables
min_lod_pixels = 5500
max_lod_pixels = -1
min_fade_extent = 0
max_fade_extent = 0



# def dms_to_dd(d, m, s):
#     dd = d + float(m)/60 + float(s)/3600
#     return dd


# def ToDMS(dd):
#         dd1 = abs(float(dd))
#         cdeg = int(dd1)
#         minsec = dd1 - cdeg
#         cmin = int(minsec * 60)
#         csec = (minsec % 60) / float(3600)
#         if dd < 0: cdeg = cdeg * -1
#         return str(cdeg,cmin,csec
# print(ToDMS(18.43348))

#npu._polygon_utils.points2convex_hull()
print("convex hull?")
print(npu.xycoords[:])
print(_polygon_utils.points2convex_hull(list(npu.xycoords[:])))


# make the polygon to display when zoomed out.
pol = kml.newpolygon(name=str(survey_title),
                     outerboundaryis=[(west_extent, north_extent),
                                    (east_extent, north_extent),
                                      (east_extent, south_extent), (west_extent, south_extent),
                                    (west_extent, north_extent)
                                      ])

                     # [(north_extent, west_extent), (north_extent, east_extent),
                     #                  (south_extent, west_extent), (south_extent, east_extent),
                     #                  (north_extent, west_extent)])

pol.style.polystyle.color = '990000ff'  # Transparent red
pol.style.polystyle.outline = 1


#kml.Region

print(simplekml.Region)
new_folder = kml.newfolder()



new_region = simplekml.Region(latlonaltbox="<north>" + str(north_extent) + "</north>" +
                                     "<south>" + str(south_extent) + "</south>" +
                                    "<east>" + str(east_extent) + "</east>" +
                                    "<west>" + str(west_extent) + "</west>",
                                    lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
				                    "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
				                    "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
				                    "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")


#ground = new_folder.newgroundoverlay(name='GroundOverlay')
new_folder.region = new_region
#ground.region = new_region


index = 0
count = npu.point_count
while index < count:
    #do the things.
    print(npu.netcdf_dataset.variables["Obsno"][index])
    obsno = npu.netcdf_dataset.variables["Obsno"][index]
    print(npu.netcdf_dataset.variables["longitude"][index])
    print(npu.netcdf_dataset.variables["latitude"][index])
    new_folder.newpoint(name=str(npu.netcdf_dataset.variables["Obsno"][index]),
        coords=[(npu.netcdf_dataset.variables["longitude"][index], npu.netcdf_dataset.variables["latitude"][index])])
    print(index)
    index = index + 1


print(survey_title)
kml.save(survey_title + ".kml")


