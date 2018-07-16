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

# Store dataset spatial extents as python variables
west_extent = npu.netcdf_dataset.getncattr('geospatial_lon_min')
east_extent = npu.netcdf_dataset.getncattr('geospatial_lon_max')
south_extent = npu.netcdf_dataset.getncattr('geospatial_lat_min')
north_extent = npu.netcdf_dataset.getncattr('geospatial_lat_max')

# set kml region constants
# Measurement in screen pixels that represents the maximum limit of the visibility range for a given Region.
MIN_LOD_PIXELS = 1200
MAX_LOD_PIXELS = -1  # -1 the default, indicates "active to infinite size."

# Distance over which the geometry fades, from fully opaque to fully transparent.
# This ramp value, expressed in screen pixels, is applied at the minimum end of the LOD (visibility) limits.
MIN_FADE_EXTENT = 200
MAX_FADE_EXTENT = 800

# point constants
POINT_ICON_STYLE_LINK = "http://maps.google.com/mapfiles/kml/shapes/placemark_square.png"

# print(west_extent)
# print("MASK")
# mask = np.zeros(shape=(npu.point_count,), dtype='bool')
# mask[-10:] = True
# print(mask)

# -----------------------bounds------------
#        self.bounds = [xmin, ymin, xmax, ymax]

def get_basic_density_measurement():
    """
    Basic calculation to get an idea of the dataset point density. This could be used to dynamically set the
    min_lod_pixels, max_lod_pixels, and point scale to get the best visual outcome. A higher density may need a smaller
    point scale or even an additional region.
    :return:
    """
    print("Point Count: {}".format(npu.point_count))
    print("Area: {}".format((east_extent - west_extent) * (north_extent - south_extent)))
    print("Density: {}".format(((east_extent - west_extent) * (north_extent - south_extent)) / npu.point_count * 1000))


def build_region(min_lod_pixels, max_lod_pixels, min_fade_extent, max_fade_extent):
    """
    Builds a kml region based on the input parameters.
    :return: the region object.
    """
    region = simplekml.Region(latlonaltbox="<north>" + str(north_extent) + "</north>" +
                                           "<south>" + str(south_extent) + "</south>" +
                                           "<east>" + str(east_extent) + "</east>" +
                                           "<west>" + str(west_extent) + "</west>",
                              lod="<minLodPixels>" + str(min_lod_pixels) + "</minLodPixels>" +
                                  "<maxLodPixels>" + str(max_lod_pixels) + "</maxLodPixels>" +
                                  "<minFadeExtent>" + str(min_fade_extent) + "</minFadeExtent>" +
                                  "<maxFadeExtent>" + str(max_fade_extent) + "</maxFadeExtent>")

    return region


def build_polygon():
    """
    Sets a containing folder for the polygon, builds a new convex hull polygon, and styles the polygon.
    This polygon is the
    :return: the folder containing the polygon, and the polygon object.
    """
    polygon_folder = kml.newfolder()

    # use _polygon_utils to find the convex hull of the netcdf. Use these points to generate a polygon.
    pol = polygon_folder.newpolygon(name=str(survey_title),
                         outerboundaryis=_polygon_utils.points2convex_hull(list(npu.xycoords[:])))

    # polygon styling
    pol.style.polystyle.color = '990000ff'  # Transparent red
    pol.style.polystyle.outline = 1

    return polygon_folder, pol


def plot_points_in_kml():
    """
    Loop through points in the netcdf file, and build and style each point as a kml place mark point.
    :return: The point folder. This folder is called to insert the corresponding region
    """
    index = 0
    count = npu.point_count
    points_folder = kml.newfolder()
    while index < count:
        # add new points with netcdf file Obsno as title and long and lat as coordinatess
        new_point = points_folder.newpoint(name=str(npu.netcdf_dataset.variables["Obsno"][index]),
        coords=[(npu.netcdf_dataset.variables["longitude"][index], npu.netcdf_dataset.variables["latitude"][index])])

        # set the point icon. Different urls can be found in point style options in google earth
        new_point.style.iconstyle.icon.href = POINT_ICON_STYLE_LINK
        new_point.style.iconstyle.scale = 0.7
        new_point.labelstyle.scale = 0  # removes the label

        index = index + 1
    assert index == count
    return points_folder

def plot_dynamic_points_in_kml():
    spatial_mask = npu.get_spatial_mask(npu.bounds)
    field_list = ['latitude', 'longitude']
    point_data_generator = npu.all_point_data_generator(field_list, spatial_mask)

    variable_attributes = next(point_data_generator)
    # print(variable_attributes.keys())
    # print(variable_attributes)

    # Use long names instead of variable names where they exist
    field_names = [variable_attributes[variable_name].get('long_name') or variable_name for variable_name in
                   variable_attributes.keys()]
    print('field_names: {}'.format(field_names))

    points_folder = kml.newfolder()
    for point_data in point_data_generator:
        print(point_data[0])
        # add new points with netcdf file Obsno as title and long and lat as coordinatess
        new_point = points_folder.newpoint(name="point", coords=[(point_data[0], point_data[1])])
        # set the point icon. Different urls can be found in point style options in google earth
        new_point.style.iconstyle.icon.href = POINT_ICON_STYLE_LINK
        new_point.style.iconstyle.scale = 0.7
        new_point.labelstyle.scale = 0  # removes the label
    print(points_folder)
    return points_folder




def build_dynamic_network_link(containing_folder):
    """
    Build a network link, set the parameters, and place in the specified folder.
    """
    net_link = containing_folder.newnetworklink(name="Network Link")
    net_link.link.href = "http://127.0.0.1:5000/query"
    net_link.link.viewrefreshmode = simplekml.ViewRefreshMode.onstop
    net_link.link.viewrefreshtime = 1
    net_link.link.refreshinterval = 2
    return net_link

def build_multiple_kmls(list_of_surveys):

    pass

def find_subset():

    lats = npu.netcdf_dataset.variables['latitude'][:]
    longs = npu.netcdf_dataset.variables['longitude'][:]
    #print(longs)
    latselect = np.logical_and(lats > 38.03, lats < 40)
    #print(latselect)
    lonselect = np.logical_and(longs > 146.6, longs < 147)
    #print(lonselect)

    #print([i for i in latselect if i is True])
    data = npu.netcdf_dataset.variables['obsno']
    #data = npu.netcdf_dataset.variables['Obsno'][0, 0, latselect, lonselect]
    #print(data)


    # lats = npu.netcdf_dataset.variables['latitude'][:]
    # lons = npu.netcdf_dataset.variables['longitude'][:]
    # latbounds = [south_extent + 0.5, north_extent - 0.5]
    # lonbounds = [west_extent + 0.5, east_extent - 0.5]  # degrees east ?
    #
    # # latitude lower and upper index
    # latli = np.argmin(np.abs(lats - latbounds[0]))
    # latui = np.argmin(np.abs(lats - latbounds[1]))
    #
    # # longitude lower and upper index
    # lonli = np.argmin(np.abs(lons - lonbounds[0]))
    # lonui = np.argmin(np.abs(lons - lonbounds[1]))
    #
    # # Air (time, latitude, longitude)
    # print(npu.netcdf_dataset.variables['Obsno'][lonli:lonui])
    # airSubset = npu.netcdf_dataset.variables['Obsno'][:, latli:latui, lonli:lonui]
    # print(airSubset)
    # print([south_extent, north_extent, west_extent, east_extent])

def main():
    #"Poorly controlled data which should be used cautiously."
    #"Data with weak gravity, position and elevation control."
    ##new_point.description = "description {0} here {1}".format(point_data[2], point_data[3])
    #get_Var = npu.get_lookup_mask(lookup_value_list=[] ,lookup_variable_name='reliab')
    #print(get_Var)
    # build the polygon, points, network link, and regions
    polygon_folder, polygon = build_polygon()
    #dataset_points_folder = plot_points_in_kml()
    #dataset_points_folder = plot_dynamic_points_in_kml()
    dataset_polygon_region = build_region(MAX_LOD_PIXELS, MIN_LOD_PIXELS + 400, MIN_FADE_EXTENT, MAX_FADE_EXTENT)
    dataset_points_region = build_region(MIN_LOD_PIXELS - 200, MAX_LOD_PIXELS, MIN_FADE_EXTENT, MAX_FADE_EXTENT)
    points_folder = kml.newfolder()
    dataset_network_link = build_dynamic_network_link(points_folder)

    # structure them correctly
    polygon_folder.region = dataset_polygon_region  # insert built polygon region into polygon folder
    #dataset_points_folder.region = dataset_points_region  # insert built point region into point folder
    dataset_network_link.region = dataset_points_region
    #find_subset()
    print("Building kml for survey: " + survey_title)
    kml.save(survey_title + ".kml")


if __name__ == '__main__':
    main()

