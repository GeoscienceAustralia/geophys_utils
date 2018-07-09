import simplekml
import netCDF4
import numpy as np
from geophys_utils import NetCDFPointUtils, get_spatial_ref_from_wkt
from geophys_utils import points2convex_hull
from geophys_utils import _polygon_utils


# Create NetCDFPointUtils object for specified netCDF dataset
netcdf_path = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/axi547/ground_gravity/point_datasets/195107.nc'
#netcdf_path = 'E:\\Temp\\gravity_point_test\\195256.nc'
#195256
#netcdf_path = "C:\\Users\\u62231\\Desktop\\grav_netcdf_4\\201780.nc"
#195105
#201780

netcdf_dataset = netCDF4.Dataset(netcdf_path)
npu = NetCDFPointUtils(netcdf_dataset)

survey_title = str(npu.netcdf_dataset.getncattr('title'))
kml = simplekml.Kml()

# Store dataset spatial extents
west_extent = npu.netcdf_dataset.getncattr('geospatial_lon_min')
east_extent = npu.netcdf_dataset.getncattr('geospatial_lon_max')
south_extent = npu.netcdf_dataset.getncattr('geospatial_lat_min')
north_extent = npu.netcdf_dataset.getncattr('geospatial_lat_max')

# set kml region constants
min_lod_pixels = 1200
max_lod_pixels = -1
min_fade_extent = 200
max_fade_extent = 800

# point constants

point_icon_style_link = "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"

print("Point Count: {}".format(npu.point_count))
print("Area: {}".format((east_extent - west_extent) * (north_extent - south_extent)))
print("Density: {}".format(((east_extent - west_extent) * (north_extent - south_extent)) / npu.point_count * 1000))





def build_polygon():

    polygon_folder = kml.newfolder()
    polygon_region = simplekml.Region(latlonaltbox="<north>" + str(north_extent) + "</north>" +
                                         "<south>" + str(south_extent) + "</south>" +
                                        "<east>" + str(east_extent) + "</east>" +
                                        "<west>" + str(west_extent) + "</west>",
                                        lod="<minLodPixels>" + str(max_lod_pixels) + "</minLodPixels>" +
                                        "<maxLodPixels>" + str(min_lod_pixels + 400) + "</maxLodPixels>" +
                                        "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
                                        "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")
    polygon_folder.region = polygon_region

    # use _polygon_utils to find the convex hull of the netcdf. Use these points to generate a polygon.
    pol = polygon_folder.newpolygon(name=str(survey_title),
                         outerboundaryis=_polygon_utils.points2convex_hull(list(npu.xycoords[:])))

    # polygon styling
    pol.style.polystyle.color = '990000ff'  # Transparent red
    pol.style.polystyle.outline = 1
    return pol


def plot_points_in_kml():

    index = 0
    count = npu.point_count
    points_folder = kml.newfolder()
    while index < count:
        # add points with Obsno as title as long and lat as coords
        new_point = points_folder.newpoint(name=str(npu.netcdf_dataset.variables["Obsno"][index]),
        coords=[(npu.netcdf_dataset.variables["longitude"][index], npu.netcdf_dataset.variables["latitude"][index])])

        # set the point icon. Different urls can be found in point style options in google earth
        new_point.style.iconstyle.icon.href = point_icon_style_link
        new_point.style.iconstyle.scale = 0.5
        new_point.labelstyle.scale = 0 # remove the label

        index = index + 1
    assert index == count
    return points_folder


def build_region(points_folder):

    points_region = simplekml.Region(latlonaltbox="<north>" + str(north_extent) + "</north>" +
                                         "<south>" + str(south_extent) + "</south>" +
                                        "<east>" + str(east_extent) + "</east>" +
                                        "<west>" + str(west_extent) + "</west>",
                                        lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
                                        "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
                                        "<minFadeExtent>" + str(0) + "</minFadeExtent>" +
                                        "<maxFadeExtent>" + str(0) + "</maxFadeExtent>")
    points_folder.region = points_region

    return points_region


def main():
    build_polygon()
    dataset_points_folder = plot_points_in_kml()
    build_region(dataset_points_folder)
    print(survey_title)
    kml.save(survey_title + ".kml")


if __name__ == '__main__':
    main()

