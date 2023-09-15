# geophys\_utils - Utilities for Accessing netCDF Encoded Geophysics Data
The geophys_utils Python library is a collection of utilities for discovering, accessing and using geophysics data via 
web services or from netCDF files. This library is maintained by Geoscience Australia with ongoing contributions from
the original Author, Alex Ip.

An overview of the netCDF encodings for grid, point and line geophysical data can be viewed at 
<https://doi.org/10.1080/22020586.2019.12073191>

Any feedback, suggestions or contributions will be gratefully accepted.

# Binder deployment
Click this button to launch an example notebook in the ARDC Binder service:

[![Binder](https://mybinder.org/badge_logo.svg)](https://binderhub.rc.nectar.org.au/v2/gh/alex-ip/geophys_utils/binder_docker_test?labpath=geophys_utils%2Fexamples%2F10_gravity_point_discovery_and_access_demo.ipynb)

Note that the BinderHub deployment uses a prebuilt base image providing a Conda environment with GDAL, netCDF4 and 
other dependencies to speed up the overall build process. The Dockerfile in the root of this repository is used to 
create the notebook container for geophys_utils in BinderHub.

# Docker deployment
There is a docker directory in the project root which contains Dockerfiles which will allow you to create fresh Docker 
images from the current source code using the command
```docker build --progress=plain -t <image_name>:<image_tag> -f docker/Dockerfile_<version> .``` 
in the project root directory. Alternatively, you can run the ```build_all_images.sh``` script to build all the images
defined.
Note that, depending on the OS, Python, and GDAL version, some builds can initially take longer due to the extensive 
dependencies for large packages such as NetCDF and SciPy. After the initial builds, the layers should be cached for 
faste subsequent builds.

Alternatively, you can pull a pre-made Ubuntu image with the latest versions using the following commands:
```bash
# Pull down the latest stable Ubuntu image from DockerHub
docker pull alexip/geophys_utils:latest
# Run unit tests
docker run --rm --name geophys_utils alexip/geophys_utils:latest python -m geophys_utils.test
# Start an interactive shell in the container
docker run -it --name geophys_utils alexip/geophys_utils:latest /bin/bash
```

## License
The content of this repository is licensed for use under the 
[Apache 2.0 License](http://www.apache.org/licenses/LICENSE-2.0). 
See the [license deed](https://github.com/GeoscienceAustralia/geophys_utils/blob/master/LICENSE) for full details.

## Contacts
**Geoscience Australia**
*Publisher*
<https://www.ga.gov.au>
[clientservices@ga.gov.au](mailto://clientservices@ga.gov.au)

**Alex Ip**
*Principal Author*
[alex@trentham.net.au](mailto://alex@trentham.net.au)
<https://orcid.org/0000-0001-8937-8904>

**Andrew Turner**
*Contributing Author*
<https://orcid.org/0000-0001-5085-8783>
